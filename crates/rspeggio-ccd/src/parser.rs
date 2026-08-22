// crates/rspeggio-ccd/src/parser.rs

use crate::component::{BondOrder, CcdAtom, CcdBond, CcdComponent};
use std::collections::HashMap;

fn load_ccd_file(path: &str) -> Result<String, std::io::Error> {
    std::fs::read_to_string(path)
}

// Parses a full CCD mmCIF file's contents into a `CcdComponent`. The atom
// loop is required — every component has atoms — so its absence or any
// malformed atom row fails the whole parse. The bond loop is optional:
// single-atom components (e.g. metal ions like Zn2+) legitimately have no
// `_chem_comp_bond` block at all, which means zero bonds, not a parse
// failure. A malformed bond *row* within a bond loop that does exist is
// still a hard failure, same as for atoms.
pub fn parse_ccd_component(contents: &str) -> Option<CcdComponent> {
    let lines: Vec<&str> = contents.lines().collect();

    let atom_block = find_loop_block(&lines, "chem_comp_atom")?;
    let atom_headers = parse_loop_headers(&atom_block.header_refs(), "chem_comp_atom");
    let atoms = build_atoms(&atom_headers, &atom_block.data_refs())?;

    let bonds = match find_loop_block(&lines, "chem_comp_bond") {
        Some(bond_block) => {
            let bond_headers = parse_loop_headers(&bond_block.header_refs(), "chem_comp_bond");
            build_bonds(&bond_headers, &bond_block.data_refs())?
        }
        None => Vec::new(),
    };

    Some(CcdComponent::new(atoms, bonds))
}

// Loads and parses a CCD mmCIF file from disk in one step. `None` covers
// both an unreadable file and a file that doesn't parse as a valid
// component -- callers needing to distinguish the two should call
// `load_ccd_file`/`parse_ccd_component` directly instead.
pub fn load_ccd_component(path: &str) -> Option<CcdComponent> {
    let contents = load_ccd_file(path).ok()?;
    parse_ccd_component(&contents)
}

// Maps each `_{category}.field_name` header line to its zero-based column
// index, so row-parsing can look fields up by name instead of position.
// The CCD schema drifts across remediation eras (extra/missing/reordered
// fields between entries), so column position alone can't be trusted.
fn parse_loop_headers(lines: &[&str], category: &str) -> HashMap<String, usize> {
    let prefix = format!("_{category}.");
    let mut headers = HashMap::new();
    for (i, line) in lines.iter().enumerate() {
        if let Some(field) = line.trim().strip_prefix(&prefix) {
            headers.insert(field.to_string(), i);
        }
    }
    headers
}

// The header lines and data lines belonging to one category's block, so a
// caller can hand `headers` straight to `parse_loop_headers` and iterate
// `data` for row parsing without re-scanning the file to find either half.
// Owned rather than borrowed: the scalar-row fallback below has to
// synthesize its single data row from several separate lines, so it can't
// stay a zero-copy slice of the original file the way a real `loop_` block
// can.
struct LoopBlock {
    headers: Vec<String>,
    data: Vec<String>,
}

impl LoopBlock {
    fn header_refs(&self) -> Vec<&str> {
        self.headers.iter().map(String::as_str).collect()
    }

    fn data_refs(&self) -> Vec<&str> {
        self.data.iter().map(String::as_str).collect()
    }
}

// Scans a CIF file's lines for the `_{category}.*` block, in either of
// mmCIF's two valid forms: a `loop_` table (many rows), or scalar
// `_{category}.field value` lines with no `loop_` at all (exactly one row —
// e.g. every monatomic ion's `_chem_comp_atom` category, which has nothing
// to bond and so also lacks any `_chem_comp_bond` block whatsoever).
// Returns `None` if the category isn't present in either form.
fn find_loop_block(lines: &[&str], category: &str) -> Option<LoopBlock> {
    let prefix = format!("_{category}.");
    find_looped_block(lines, &prefix).or_else(|| find_scalar_block(lines, &prefix))
}

fn find_looped_block(lines: &[&str], prefix: &str) -> Option<LoopBlock> {
    let mut i = 0;
    while i < lines.len() {
        if lines[i].trim() == "loop_" {
            let header_start = i + 1;
            if header_start < lines.len() && lines[header_start].trim().starts_with(prefix) {
                let mut header_end = header_start;
                while header_end < lines.len() && lines[header_end].trim().starts_with(prefix) {
                    header_end += 1;
                }
                let mut data_end = header_end;
                while data_end < lines.len() {
                    let t = lines[data_end].trim();
                    if t.is_empty() || t.starts_with('#') || t == "loop_" || t.starts_with('_') {
                        break;
                    }
                    data_end += 1;
                }
                return Some(LoopBlock {
                    headers: lines[header_start..header_end]
                        .iter()
                        .map(|s| s.to_string())
                        .collect(),
                    data: lines[header_end..data_end]
                        .iter()
                        .map(|s| s.to_string())
                        .collect(),
                });
            }
        }
        i += 1;
    }
    None
}

// Collects `_{category}.field value` lines and synthesizes them into one
// data row, in the order the fields appeared, so the rest of the pipeline
// (header map + tokenize-a-row) works identically regardless of which form
// the file used.
fn find_scalar_block(lines: &[&str], prefix: &str) -> Option<LoopBlock> {
    let mut headers = Vec::new();
    let mut values = Vec::new();
    for line in lines {
        let Some(rest) = line.trim().strip_prefix(prefix) else {
            continue;
        };
        let mut parts = rest.splitn(2, char::is_whitespace);
        let field = parts.next()?;
        let value = parts.next().unwrap_or("").trim();
        headers.push(format!("{prefix}{field}"));
        values.push(value.to_string());
    }
    if headers.is_empty() {
        return None;
    }
    Some(LoopBlock {
        headers,
        data: vec![values.join(" ")],
    })
}

// Builds `CcdAtom`s from a `_chem_comp_atom` loop's header map and data
// rows. Returns `None` if a required field is missing from the header map,
// or if any data row is short a token for one of those fields — CCD schema
// drift is expected across categories/eras, but a malformed atom row within
// a category we *did* find is a real problem, not something to paper over.
fn build_atoms(headers: &HashMap<String, usize>, data: &[&str]) -> Option<Vec<CcdAtom>> {
    let atom_id_idx = *headers.get("atom_id")?;
    let element_idx = *headers.get("type_symbol")?;
    let aromatic_idx = *headers.get("pdbx_aromatic_flag")?;
    let leaving_idx = *headers.get("pdbx_leaving_atom_flag")?;

    data.iter()
        .map(|line| {
            let tokens = tokenize_cif_row(line);
            let atom_id = tokens.get(atom_id_idx)?.clone();
            let element = tokens.get(element_idx)?.clone();
            let aromatic = tokens.get(aromatic_idx)?.as_str() == "Y";
            let leaving = tokens.get(leaving_idx)?.as_str() == "Y";
            Some(CcdAtom::new(atom_id, element, aromatic, leaving))
        })
        .collect()
}

// Builds `CcdBond`s from a `_chem_comp_bond` loop's header map and data
// rows. Same failure posture as `build_atoms`: a missing required header,
// a short row, or an unrecognized `value_order` token all yield `None`
// rather than a guessed bond.
fn build_bonds(headers: &HashMap<String, usize>, data: &[&str]) -> Option<Vec<CcdBond>> {
    let atom_1_idx = *headers.get("atom_id_1")?;
    let atom_2_idx = *headers.get("atom_id_2")?;
    let order_idx = *headers.get("value_order")?;
    let aromatic_idx = *headers.get("pdbx_aromatic_flag")?;

    data.iter()
        .map(|line| {
            let tokens = tokenize_cif_row(line);
            let atom_id_1 = tokens.get(atom_1_idx)?.clone();
            let atom_id_2 = tokens.get(atom_2_idx)?.clone();
            let order = BondOrder::from_ccd_str(tokens.get(order_idx)?)?;
            let aromatic = tokens.get(aromatic_idx)?.as_str() == "Y";
            Some(CcdBond::new(atom_id_1, atom_id_2, order, aromatic))
        })
        .collect()
}

// TODO: in the future this could be zero-copy borrowed return <Vec &str> but for now will copy
// TODO: revisit this as an iterator-based walk (chars.iter().peekable(), .next()/.peek())
// instead of manual index bookkeeping, once quote-handling (step 3) is done and this is
// fully working end to end. Not a correctness issue, just a more idiomatic-Rust refactor
// to come back to deliberately, not mid-implementation.
fn tokenize_cif_row(line: &str) -> Vec<String> {
    let chars: Vec<char> = line.chars().collect();
    let mut tokens: Vec<String> = Vec::new();
    let mut i = 0;
    let mut start: usize = 0;
    while i < chars.len() {
        // 1. skip any whitespace characters, advancing i
        // 2. if we've hit the end of the line, stop
        // 3. otherwise, read characters until the next whitespace (or end of line)
        // 4. collect that range of chars into a String, push onto tokens
        if chars[i] == '"' || chars[i] == '\'' {
            let quote_char = chars[i];
            let mut j = i + 1;
            while j < chars.len() {
                if chars[j] == quote_char && (j + 1 >= chars.len() || chars[j + 1].is_whitespace())
                {
                    break;
                }
                j += 1;
            }
            let token: String = chars[i + 1..j].iter().collect();
            tokens.push(token);
            i = j + 1;
            start = i;
        } else if chars[i].is_whitespace() {
            if start < i {
                let token: String = chars[start..i].iter().collect(); // not inclusive of i
                tokens.push(token);
            }
            start = i + 1;
            i += 1;
        } else {
            i += 1;
        }
    }
    if start < chars.len() {
        let final_token: String = chars[start..chars.len()].iter().collect();
        tokens.push(final_token);
    }
    tokens
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn loads_a_real_ccd_file() {
        let result = load_ccd_file("tests/fixtures/ADP_ideal.cif");
        assert!(result.is_ok());

        let contents = result.expect("fixture file should exist and be readable");
        assert!(contents.contains("ADP"));
    }

    #[test]
    fn empty_line_yields_no_tokens() {
        assert_eq!(tokenize_cif_row(""), Vec::<String>::new());
    }

    #[test]
    fn splits_unquoted_fields_on_whitespace() {
        let tokens = tokenize_cif_row("ALA N N N 0 1 N N N N N N");
        assert_eq!(
            tokens,
            vec!["ALA", "N", "N", "N", "0", "1", "N", "N", "N", "N", "N", "N"]
        );
    }

    #[test]
    fn maps_atom_header_fields_to_column_indices() {
        let lines = [
            "_chem_comp_atom.comp_id ",
            "_chem_comp_atom.atom_id ",
            "_chem_comp_atom.alt_atom_id ",
            "_chem_comp_atom.type_symbol ",
            "_chem_comp_atom.charge ",
            "_chem_comp_atom.pdbx_align ",
            "_chem_comp_atom.pdbx_aromatic_flag ",
            "_chem_comp_atom.pdbx_leaving_atom_flag ",
            "_chem_comp_atom.pdbx_stereo_config ",
            "_chem_comp_atom.pdbx_backbone_atom_flag ",
            "_chem_comp_atom.pdbx_n_terminal_atom_flag ",
            "_chem_comp_atom.pdbx_c_terminal_atom_flag ",
            "_chem_comp_atom.model_Cartn_x ",
            "_chem_comp_atom.model_Cartn_y ",
            "_chem_comp_atom.model_Cartn_z ",
            "_chem_comp_atom.pdbx_model_Cartn_x_ideal ",
            "_chem_comp_atom.pdbx_model_Cartn_y_ideal ",
            "_chem_comp_atom.pdbx_model_Cartn_z_ideal ",
            "_chem_comp_atom.pdbx_component_atom_id ",
            "_chem_comp_atom.pdbx_component_comp_id ",
            "_chem_comp_atom.pdbx_ordinal ",
        ];

        let headers = parse_loop_headers(&lines, "chem_comp_atom");

        assert_eq!(headers.len(), lines.len());
        assert_eq!(headers.get("comp_id"), Some(&0));
        assert_eq!(headers.get("atom_id"), Some(&1));
        assert_eq!(headers.get("type_symbol"), Some(&3));
        assert_eq!(headers.get("pdbx_aromatic_flag"), Some(&6));
        assert_eq!(headers.get("pdbx_ordinal"), Some(&20));
    }

    #[test]
    fn header_parsing_is_reusable_across_loop_categories() {
        let lines = [
            "_chem_comp_bond.comp_id ",
            "_chem_comp_bond.atom_id_1 ",
            "_chem_comp_bond.atom_id_2 ",
            "_chem_comp_bond.value_order ",
            "_chem_comp_bond.pdbx_aromatic_flag ",
            "_chem_comp_bond.pdbx_stereo_config ",
            "_chem_comp_bond.pdbx_ordinal ",
        ];

        let headers = parse_loop_headers(&lines, "chem_comp_bond");

        assert_eq!(headers.len(), 7);
        assert_eq!(headers.get("atom_id_1"), Some(&1));
        assert_eq!(headers.get("atom_id_2"), Some(&2));
        assert_eq!(headers.get("pdbx_ordinal"), Some(&6));
        // fields belonging to a different category are never present
        assert_eq!(headers.get("type_symbol"), None);
    }

    #[test]
    fn empty_header_lines_yield_empty_map() {
        let lines: [&str; 0] = [];
        assert!(parse_loop_headers(&lines, "chem_comp_atom").is_empty());
    }

    #[test]
    fn finds_the_atom_loop_block_in_a_real_ccd_file() {
        let contents = load_ccd_file("tests/fixtures/ADP_ideal.cif").expect("fixture should load");
        let lines: Vec<&str> = contents.lines().collect();

        let block = find_loop_block(&lines, "chem_comp_atom").expect("atom loop should be found");

        assert_eq!(block.headers.len(), 21);
        assert_eq!(block.headers[0].trim(), "_chem_comp_atom.comp_id");
        assert_eq!(block.headers[20].trim(), "_chem_comp_atom.pdbx_ordinal");

        assert_eq!(block.data.len(), 42);
        assert!(block.data[0].trim_start().starts_with("ADP PB"));
        assert!(block
            .data
            .last()
            .unwrap()
            .trim_start()
            .starts_with("ADP H2"));
    }

    #[test]
    fn finds_the_bond_loop_block_in_a_real_ccd_file() {
        let contents = load_ccd_file("tests/fixtures/ADP_ideal.cif").expect("fixture should load");
        let lines: Vec<&str> = contents.lines().collect();

        let block = find_loop_block(&lines, "chem_comp_bond").expect("bond loop should be found");

        assert_eq!(block.headers.len(), 7);
        assert_eq!(block.data.len(), 44);
        // atom and bond data slices must not overlap
        let atom_block = find_loop_block(&lines, "chem_comp_atom").unwrap();
        assert!(atom_block.data.last().unwrap() != block.data.first().unwrap());
    }

    #[test]
    fn missing_category_returns_none() {
        let contents = load_ccd_file("tests/fixtures/ADP_ideal.cif").expect("fixture should load");
        let lines: Vec<&str> = contents.lines().collect();

        assert!(find_loop_block(&lines, "chem_comp_nonexistent").is_none());
    }

    #[test]
    fn builds_atoms_from_the_real_ccd_fixture() {
        let contents = load_ccd_file("tests/fixtures/ADP_ideal.cif").expect("fixture should load");
        let lines: Vec<&str> = contents.lines().collect();
        let block = find_loop_block(&lines, "chem_comp_atom").expect("atom loop should be found");
        let headers = parse_loop_headers(&block.header_refs(), "chem_comp_atom");

        let atoms = build_atoms(&headers, &block.data_refs()).expect("atoms should build");

        assert_eq!(atoms.len(), 42);
        assert_eq!(
            atoms[0],
            CcdAtom::new("PB".to_string(), "P".to_string(), false, false)
        );
        let n9 = atoms
            .iter()
            .find(|a| a.atom_id == "N9")
            .expect("N9 atom should be present");
        assert_eq!(n9.element, "N");
        assert!(n9.aromatic, "N9 is part of the purine ring, flagged Y");
        assert!(
            !n9.leaving,
            "N9 is a core ring atom, not removed on polymerization"
        );
    }

    #[test]
    fn build_atoms_fails_when_a_required_header_is_missing() {
        let mut headers = HashMap::new();
        headers.insert("atom_id".to_string(), 0);
        // "type_symbol" deliberately absent
        headers.insert("pdbx_aromatic_flag".to_string(), 2);
        headers.insert("pdbx_leaving_atom_flag".to_string(), 3);

        assert!(build_atoms(&headers, &["PB P N N"]).is_none());
    }

    #[test]
    fn build_atoms_fails_on_a_short_data_row() {
        let mut headers = HashMap::new();
        headers.insert("atom_id".to_string(), 0);
        headers.insert("type_symbol".to_string(), 1);
        headers.insert("pdbx_aromatic_flag".to_string(), 2);
        headers.insert("pdbx_leaving_atom_flag".to_string(), 3);

        // only two tokens, but pdbx_aromatic_flag is expected at index 2
        assert!(build_atoms(&headers, &["PB P"]).is_none());
    }

    #[test]
    fn builds_bonds_from_the_real_ccd_fixture() {
        let contents = load_ccd_file("tests/fixtures/ADP_ideal.cif").expect("fixture should load");
        let lines: Vec<&str> = contents.lines().collect();
        let block = find_loop_block(&lines, "chem_comp_bond").expect("bond loop should be found");
        let headers = parse_loop_headers(&block.header_refs(), "chem_comp_bond");

        let bonds = build_bonds(&headers, &block.data_refs()).expect("bonds should build");

        assert_eq!(bonds.len(), 44);
        assert_eq!(
            bonds[0],
            CcdBond::new(
                "PB".to_string(),
                "O1B".to_string(),
                BondOrder::Double,
                false
            )
        );
        assert_eq!(
            bonds[1],
            CcdBond::new(
                "PB".to_string(),
                "O2B".to_string(),
                BondOrder::Single,
                false
            )
        );
    }

    #[test]
    fn parses_a_full_component_from_a_real_ccd_file() {
        let contents = load_ccd_file("tests/fixtures/ADP_ideal.cif").expect("fixture should load");

        let component = parse_ccd_component(&contents).expect("ADP should parse");

        assert_eq!(component.atoms().len(), 42);
        assert_eq!(component.bonds().len(), 44);
    }

    #[test]
    fn finds_a_scalar_form_atom_block_with_no_loop_keyword() {
        let contents = load_ccd_file("tests/fixtures/ZN_ideal.cif").expect("fixture should load");
        let lines: Vec<&str> = contents.lines().collect();

        let block = find_loop_block(&lines, "chem_comp_atom")
            .expect("scalar-form atom block should be found");

        assert_eq!(
            block.data.len(),
            1,
            "scalar form synthesizes exactly one row"
        );
        let headers = parse_loop_headers(&block.header_refs(), "chem_comp_atom");
        let atoms = build_atoms(&headers, &block.data_refs()).expect("atom should build");

        assert_eq!(atoms.len(), 1);
        assert_eq!(atoms[0].atom_id, "ZN");
        assert_eq!(atoms[0].element, "ZN");
        assert!(!atoms[0].aromatic);
        assert!(
            !atoms[0].leaving,
            "a bare ion has nothing to leave on polymerization"
        );
    }

    #[test]
    fn loads_a_component_straight_from_a_file_path() {
        let component =
            load_ccd_component("tests/fixtures/ADP_ideal.cif").expect("ADP should load");

        assert_eq!(component.atoms().len(), 42);
        assert_eq!(component.bonds().len(), 44);
    }

    #[test]
    fn load_ccd_component_returns_none_for_a_missing_file() {
        assert!(load_ccd_component("tests/fixtures/does_not_exist.cif").is_none());
    }

    #[test]
    fn parses_a_bondless_single_atom_component() {
        let contents = load_ccd_file("tests/fixtures/ZN_ideal.cif").expect("fixture should load");

        let component = parse_ccd_component(&contents).expect("ZN should parse");

        assert_eq!(component.atoms().len(), 1);
        assert_eq!(component.atoms()[0].atom_id, "ZN");
        // no _chem_comp_bond loop at all for a bare ion -- zero bonds, not a failure
        assert!(component.bonds().is_empty());
    }

    #[test]
    fn parses_a_component_whose_separator_line_is_hash_space_hash() {
        // PO4's real RCSB export uses "#   #" as its loop separator (a
        // fixed-width export quirk) instead of ADP's plain "# " -- schema
        // drift that once made `find_loop_block` swallow the separator as
        // a malformed data row and fail the whole parse.
        let component = load_ccd_component("data/common/PO4.cif").expect("PO4 should parse");

        assert_eq!(component.atoms().len(), 5, "P + O1..O4");
        assert_eq!(component.bonds().len(), 4, "P-O1..P-O4");
    }

    #[test]
    fn builds_a_real_triple_bond_from_the_cyanide_fixture() {
        let contents = load_ccd_file("tests/fixtures/CN_ideal.cif").expect("fixture should load");
        let lines: Vec<&str> = contents.lines().collect();
        let block = find_loop_block(&lines, "chem_comp_bond").expect("bond loop should be found");
        let headers = parse_loop_headers(&block.header_refs(), "chem_comp_bond");

        let bonds = build_bonds(&headers, &block.data_refs()).expect("bonds should build");

        assert_eq!(bonds.len(), 2);
        assert_eq!(
            bonds[0],
            CcdBond::new("C1".to_string(), "N1".to_string(), BondOrder::Triple, false)
        );
        assert_eq!(
            bonds[1],
            CcdBond::new("C1".to_string(), "H1".to_string(), BondOrder::Single, false)
        );
    }

    #[test]
    fn build_bonds_fails_when_a_required_header_is_missing() {
        let mut headers = HashMap::new();
        headers.insert("atom_id_1".to_string(), 0);
        headers.insert("atom_id_2".to_string(), 1);
        // "value_order" deliberately absent
        headers.insert("pdbx_aromatic_flag".to_string(), 3);

        assert!(build_bonds(&headers, &["PB O1B DOUB N"]).is_none());
    }

    #[test]
    fn build_bonds_fails_on_an_unrecognized_value_order() {
        let mut headers = HashMap::new();
        headers.insert("atom_id_1".to_string(), 0);
        headers.insert("atom_id_2".to_string(), 1);
        headers.insert("value_order".to_string(), 2);
        headers.insert("pdbx_aromatic_flag".to_string(), 3);

        assert!(build_bonds(&headers, &["PB O1B QUAD N"]).is_none());
    }

    #[test]
    fn bond_order_recognizes_all_four_ccd_values() {
        assert_eq!(BondOrder::from_ccd_str("SING"), Some(BondOrder::Single));
        assert_eq!(BondOrder::from_ccd_str("DOUB"), Some(BondOrder::Double));
        assert_eq!(BondOrder::from_ccd_str("TRIP"), Some(BondOrder::Triple));
        assert_eq!(BondOrder::from_ccd_str("AROM"), Some(BondOrder::Aromatic));
        assert_eq!(BondOrder::from_ccd_str("QUAD"), None);
    }

    #[test]
    fn handles_quoted_fields_with_apostrophes() {
        let tokens = tokenize_cif_row("TPP \"N1'\"  \"N1'\"  N 0 1 Y N N");
        assert_eq!(
            tokens,
            vec!["TPP", "N1'", "N1'", "N", "0", "1", "Y", "N", "N"]
        );
    }
}
