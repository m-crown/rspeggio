// crates/rspeggio-ccd/src/parser.rs
// TODO: remove once load_ccd_file/tokenize_cif_row are called from real (non-test) code,
// and CcdAtom/CcdComponent are actually constructed here.
#![allow(dead_code)]
#![allow(unused_imports)]

use crate::component::{CcdAtom, CcdComponent};
use std::collections::HashMap;

fn load_ccd_file(path: &str) -> Result<String, std::io::Error> {
    std::fs::read_to_string(path)
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

// The header lines and data lines belonging to one `loop_` block, so a
// caller can hand `headers` straight to `parse_loop_headers` and iterate
// `data` for row parsing without re-scanning the file to find either half.
struct LoopBlock<'a> {
    headers: &'a [&'a str],
    data: &'a [&'a str],
}

// Scans a CIF file's lines for the `loop_` block whose headers start with
// `_{category}.`, returning `None` if that category isn't present at all
// (some CCD entries omit e.g. a bond loop).
fn find_loop_block<'a>(lines: &'a [&'a str], category: &str) -> Option<LoopBlock<'a>> {
    let prefix = format!("_{category}.");
    let mut i = 0;
    while i < lines.len() {
        if lines[i].trim() == "loop_" {
            let header_start = i + 1;
            if header_start < lines.len() && lines[header_start].trim().starts_with(&prefix) {
                let mut header_end = header_start;
                while header_end < lines.len() && lines[header_end].trim().starts_with(&prefix) {
                    header_end += 1;
                }
                let mut data_end = header_end;
                while data_end < lines.len() {
                    let t = lines[data_end].trim();
                    if t.is_empty() || t == "#" || t == "loop_" || t.starts_with('_') {
                        break;
                    }
                    data_end += 1;
                }
                return Some(LoopBlock {
                    headers: &lines[header_start..header_end],
                    data: &lines[header_end..data_end],
                });
            }
        }
        i += 1;
    }
    None
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

    data.iter()
        .map(|line| {
            let tokens = tokenize_cif_row(line);
            let atom_id = tokens.get(atom_id_idx)?.clone();
            let element = tokens.get(element_idx)?.clone();
            let aromatic = tokens.get(aromatic_idx)?.as_str() == "Y";
            Some(CcdAtom::new(atom_id, element, aromatic))
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
        let headers = parse_loop_headers(block.headers, "chem_comp_atom");

        let atoms = build_atoms(&headers, block.data).expect("atoms should build");

        assert_eq!(atoms.len(), 42);
        assert_eq!(
            atoms[0],
            CcdAtom::new("PB".to_string(), "P".to_string(), false)
        );
        let n9 = atoms
            .iter()
            .find(|a| a.atom_id == "N9")
            .expect("N9 atom should be present");
        assert_eq!(n9.element, "N");
        assert!(n9.aromatic, "N9 is part of the purine ring, flagged Y");
    }

    #[test]
    fn build_atoms_fails_when_a_required_header_is_missing() {
        let mut headers = HashMap::new();
        headers.insert("atom_id".to_string(), 0);
        // "type_symbol" deliberately absent
        headers.insert("pdbx_aromatic_flag".to_string(), 2);

        assert!(build_atoms(&headers, &["PB P N"]).is_none());
    }

    #[test]
    fn build_atoms_fails_on_a_short_data_row() {
        let mut headers = HashMap::new();
        headers.insert("atom_id".to_string(), 0);
        headers.insert("type_symbol".to_string(), 1);
        headers.insert("pdbx_aromatic_flag".to_string(), 2);

        // only two tokens, but pdbx_aromatic_flag is expected at index 2
        assert!(build_atoms(&headers, &["PB P"]).is_none());
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
