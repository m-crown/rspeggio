// crates/rspeggio-core/src/selection.rs
//
// Real pdbe-arpeggio's selection mechanism (`utils.selection_parser`,
// `Arpeggio._make_selection`): the CLI's `-s` argument, which says which
// residues/atoms are "the ligand" (or whatever's being analysed). Every
// `interacting_entities` classification in `export.rs`
// (INTER/INTRA_SELECTION/INTRA_NON_SELECTION/SELECTION_WATER/...) depends
// entirely on it -- and so does *which contacts get computed at all*: real
// arpeggio rebuilds its neighbor search over only `selection_plus` (the
// selection expanded to its binding site), so an atom pair where neither
// side is anywhere near the selection is never even considered.
//
// Two real syntaxes are ported: `RESNAME:<comp_id>` (every atom in every
// residue with that exact comp_id) and `/<chain>/<resnum>[<inscode>]/<atom>`
// (each field optional, filters progressively, real fields separated by
// `/`). `LIGANDS` (real arpeggio's heuristic auto-detector for "the" bound
// small molecule) is NOT implemented -- a real, deliberate gap, not
// silently dropped: none of this project's golden fixtures needed it.
// Their real selections were worked out by hand-checking
// `tests/fixtures/golden/*.json`'s own `interacting_entities` values (not
// guessed) and all turn out to be plain `RESNAME:<id>` selections -- ZN
// for 1CA2, HEM for 1MBO, FMN for 1FLV, FES for 4FXC.
//
// No selection at all (`SelectionContext::whole_structure`) reproduces
// real arpeggio's own documented default when `-s` is omitted entirely:
// the whole structure becomes "the selection"
// (`_make_selection`: `selection = entity if not selections else ...`).
// 1UBQ's own golden fixture was generated exactly this way (confirmed: its
// `interacting_entities` values are only INTRA_SELECTION/SELECTION_WATER/
// WATER_WATER, never INTER/INTRA_NON_SELECTION -- the same collapse this
// module's `whole_structure` derives from the general rules below, not a
// separately hand-maintained special case).

use pdbtbx::{
    Atom, AtomConformerResidueChainModel, ContainsAtomConformer, ContainsAtomConformerResidue,
    ContainsAtomConformerResidueChain, Residue, PDB,
};
use std::collections::HashSet;

// Real pdbe-arpeggio's binding-site expansion radius
// (`_make_selection`'s `self.ns.search_all(6.0)`) -- a separate real
// constant from `config::CONTACT_TYPES_MAX_DIST` (4.5A), specific to
// selection expansion, not contact detection.
pub const BINDING_SITE_EXPANSION_DISTANCE: f64 = 6.0;

fn is_digits(s: &str) -> bool {
    !s.is_empty() && s.chars().all(|c| c.is_ascii_digit())
}

#[derive(Debug, Clone)]
enum SelectionSpec {
    Resname(String),
    Path {
        chain: Option<String>,
        residue_number: Option<isize>,
        insertion_code: String,
        atom_name: Option<String>,
    },
}

// Ported from real `utils.selection_parser`'s per-token grammar. Same
// failure posture as the rest of this project's parsing code: a malformed
// or unrecognized selection string is a hard error, not silently ignored.
fn parse_spec(raw: &str) -> Result<SelectionSpec, String> {
    let spec = raw.trim();

    if let Some(rest) = spec.strip_prefix("RESNAME:") {
        let comp_id = rest.trim();
        if comp_id.is_empty() || comp_id.len() > 3 {
            return Err(format!(
                "invalid selection {raw:?}: RESNAME id must be 1-3 characters"
            ));
        }
        return Ok(SelectionSpec::Resname(comp_id.to_string()));
    }

    if let Some(rest) = spec.strip_prefix('/') {
        let fields: Vec<&str> = rest.split('/').collect();
        if fields.len() != 3 {
            return Err(format!(
                "invalid selection {raw:?}: expected /<chain>/<resnum>[<inscode>]/<atom>"
            ));
        }

        let chain = (!fields[0].is_empty()).then(|| fields[0].to_string());

        let mut residue_number = None;
        let mut insertion_code = " ".to_string();
        if !fields[1].is_empty() {
            if is_digits(fields[1]) {
                residue_number = fields[1].parse::<isize>().ok();
            } else if fields[1].chars().all(|c| c.is_ascii_alphanumeric()) {
                let (number_part, code_part) = fields[1].split_at(fields[1].len() - 1);
                let code_is_alpha = code_part
                    .chars()
                    .next()
                    .is_some_and(|c| c.is_ascii_alphabetic());
                if code_is_alpha && is_digits(number_part) {
                    residue_number = number_part.parse::<isize>().ok();
                    insertion_code = code_part.to_string();
                } else {
                    return Err(format!(
                        "invalid selection {raw:?}: bad residue number/inscode"
                    ));
                }
            } else {
                return Err(format!(
                    "invalid selection {raw:?}: bad residue number/inscode"
                ));
            }
            if residue_number.is_none() {
                return Err(format!("invalid selection {raw:?}: bad residue number"));
            }
        }

        let atom_name = (!fields[2].is_empty()).then(|| fields[2].to_string());

        return Ok(SelectionSpec::Path {
            chain,
            residue_number,
            insertion_code,
            atom_name,
        });
    }

    Err(format!(
        "invalid selection {raw:?}: expected RESNAME:<id> or /<chain>/<resnum>[<inscode>]/<atom>"
    ))
}

fn matches_spec(spec: &SelectionSpec, hierarchy: &AtomConformerResidueChainModel) -> bool {
    match spec {
        SelectionSpec::Resname(comp_id) => hierarchy.residue().name() == Some(comp_id.as_str()),
        SelectionSpec::Path {
            chain,
            residue_number,
            insertion_code,
            atom_name,
        } => {
            if let Some(chain_id) = chain {
                if hierarchy.chain().id() != chain_id {
                    return false;
                }
            }
            if let Some(number) = residue_number {
                let (real_number, real_ins) = hierarchy.residue().id();
                if real_number != *number || real_ins.unwrap_or(" ") != insertion_code {
                    return false;
                }
            }
            if let Some(name) = atom_name {
                if hierarchy.atom().name() != name {
                    return false;
                }
            }
            true
        }
    }
}

// A real atom's identity, usable as a hash key across thread boundaries.
// Real atom identity is still by pointer into the live `PDB` (same
// technique the rest of this crate already uses for residue/ring dedup,
// e.g. `export.rs`'s `residue_instances`) -- but a bare `*const Atom`
// itself is neither `Send` nor `Sync`, which M8's rayon parallelism needs
// to move these across threads. The pointer's own address, as a `usize`,
// carries the same identity and is a plain, thread-safe integer.
type AtomKey = usize;

fn atom_key(atom: &Atom) -> AtomKey {
    atom as *const Atom as usize
}

// A real, resolved set of atoms -- either an explicit user selection or
// its binding-site expansion.
pub struct Selection {
    atoms: HashSet<AtomKey>,
}

impl Selection {
    fn from_atom_keys(atoms: HashSet<AtomKey>) -> Self {
        Self { atoms }
    }

    pub fn whole_structure(pdb: &PDB) -> Self {
        Self::from_atom_keys(
            pdb.atoms_with_hierarchy()
                .map(|h| atom_key(h.atom()))
                .collect(),
        )
    }

    // `Err` if no real atom matches any spec -- real `selection_parser`
    // raises `SelectionError('entity not found')` for exactly this case.
    pub fn parse(pdb: &PDB, specs: &[String]) -> Result<Self, String> {
        let parsed: Vec<SelectionSpec> = specs
            .iter()
            .map(|s| parse_spec(s))
            .collect::<Result<_, _>>()?;

        let mut atoms = HashSet::new();
        for hierarchy in pdb.atoms_with_hierarchy() {
            if parsed.iter().any(|spec| matches_spec(spec, &hierarchy)) {
                atoms.insert(atom_key(hierarchy.atom()));
            }
        }
        if atoms.is_empty() {
            return Err("selection matched no real atoms".to_string());
        }
        Ok(Self::from_atom_keys(atoms))
    }

    pub fn contains(&self, atom: &Atom) -> bool {
        self.atoms.contains(&atom_key(atom))
    }

    // Real ring-ring/atom-plane/group-group/group-plane membership checks
    // are all at *residue* granularity (`ring['residue'] in
    // selection_residues`, where `selection_residues` is built from the
    // owning residues of every selected atom) -- a residue counts as "in"
    // if *any* of its real resolved atoms is.
    pub fn contains_any_atom_of(&self, residue: &Residue) -> bool {
        residue.atoms().any(|atom| self.contains(atom))
    }
}

// The two real atom sets every contact classification needs: the raw
// selection, and its expansion to the binding site (any real atom within
// `BINDING_SITE_EXPANSION_DISTANCE` of any selected atom). Real arpeggio
// only ever computes contacts among `selection_plus` atoms at all --
// `export.rs`'s contact-collection functions filter on this before
// classifying anything.
pub struct SelectionContext {
    pub selection: Selection,
    pub selection_plus: Selection,
}

impl SelectionContext {
    pub fn whole_structure(pdb: &PDB) -> Self {
        Self {
            selection: Selection::whole_structure(pdb),
            selection_plus: Selection::whole_structure(pdb),
        }
    }

    pub fn from_specs(pdb: &PDB, specs: &[String]) -> Result<Self, String> {
        let selection = Selection::parse(pdb, specs)?;
        let selection_plus = expand_to_binding_site(pdb, &selection);
        Ok(Self {
            selection,
            selection_plus,
        })
    }
}

fn expand_to_binding_site(pdb: &PDB, selection: &Selection) -> Selection {
    let tree = pdb.create_hierarchy_rtree();
    let cutoff_squared = BINDING_SITE_EXPANSION_DISTANCE * BINDING_SITE_EXPANSION_DISTANCE;

    let mut plus: HashSet<AtomKey> = selection.atoms.clone();
    for hierarchy in tree.iter() {
        if !selection.contains(hierarchy.atom()) {
            continue;
        }
        for neighbor in tree.locate_within_distance(hierarchy.atom().pos(), cutoff_squared) {
            plus.insert(atom_key(neighbor.atom()));
        }
    }
    Selection::from_atom_keys(plus)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resname_selection_matches_real_zinc_atoms_in_1ca2() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1CA2.cif").expect("1CA2 should load");

        let selection =
            Selection::parse(&pdb, &["RESNAME:ZN".to_string()]).expect("ZN should be found");

        let zn = pdb
            .atoms_with_hierarchy()
            .find(|h| h.residue().name() == Some("ZN"))
            .expect("1CA2 has a real ZN atom");
        assert!(selection.contains(zn.atom()));

        let non_zn = pdb
            .atoms_with_hierarchy()
            .find(|h| h.residue().name() != Some("ZN"))
            .expect("1CA2 has non-ZN atoms too");
        assert!(!selection.contains(non_zn.atom()));
    }

    #[test]
    fn a_resname_longer_than_three_characters_is_rejected() {
        assert!(parse_spec("RESNAME:GLYX").is_err());
    }

    #[test]
    fn path_selection_matches_a_real_specific_atom() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");

        let selection =
            Selection::parse(&pdb, &["/A/1/CA".to_string()]).expect("MET1 CA should be found");

        let met1_ca = pdb
            .atoms_with_hierarchy()
            .find(|h| h.residue().id().0 == 1 && h.atom().name() == "CA")
            .expect("MET1 CA should be present");
        assert!(selection.contains(met1_ca.atom()));

        let met1_n = pdb
            .atoms_with_hierarchy()
            .find(|h| h.residue().id().0 == 1 && h.atom().name() == "N")
            .expect("MET1 N should be present");
        assert!(
            !selection.contains(met1_n.atom()),
            "the atom-name field should narrow the selection to just CA"
        );
    }

    #[test]
    fn path_selection_with_omitted_fields_selects_a_whole_residue() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");

        let selection =
            Selection::parse(&pdb, &["/A/1/".to_string()]).expect("MET1 should be found");

        for atom_name in ["N", "CA", "C", "O"] {
            let atom = pdb
                .atoms_with_hierarchy()
                .find(|h| h.residue().id().0 == 1 && h.atom().name() == atom_name)
                .unwrap_or_else(|| panic!("MET1 {atom_name} should be present"));
            assert!(selection.contains(atom.atom()));
        }
    }

    #[test]
    fn an_unmatched_selection_is_an_error() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");
        assert!(Selection::parse(&pdb, &["RESNAME:ZZZ".to_string()]).is_err());
    }

    #[test]
    fn a_malformed_selection_string_is_an_error() {
        assert!(parse_spec("nonsense").is_err());
        assert!(parse_spec("/A/B/C/D").is_err());
    }

    #[test]
    fn binding_site_expansion_includes_real_nearby_atoms_but_not_far_ones() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1CA2.cif").expect("1CA2 should load");

        let context =
            SelectionContext::from_specs(&pdb, &["RESNAME:ZN".to_string()]).expect("should build");

        // NE2 specifically -- the real side-chain atom that directly
        // coordinates the zinc (~2.0A away, confirmed via an ad-hoc
        // distance check before writing this test). A different atom in
        // the same residue (e.g. the backbone N) isn't guaranteed to be
        // within 6A just because its side chain's tip is.
        let his94_ne2 = pdb
            .atoms_with_hierarchy()
            .find(|h| {
                h.residue().name() == Some("HIS")
                    && h.residue().id().0 == 94
                    && h.atom().name() == "NE2"
            })
            .expect("1CA2 has a real HIS94 NE2 (a real zinc ligand)");
        assert!(
            context.selection_plus.contains(his94_ne2.atom()),
            "HIS94's real zinc-coordinating NE2 should be within the 6A binding site expansion"
        );

        let total_atoms = pdb.atoms_with_hierarchy().count();
        let plus_count = pdb
            .atoms_with_hierarchy()
            .filter(|h| context.selection_plus.contains(h.atom()))
            .count();
        assert!(
            plus_count < total_atoms,
            "the binding site expansion around one real zinc atom should be a strict subset \
             of 1CA2's ~2000 real atoms, not the whole structure"
        );
        assert!(
            plus_count > 1,
            "expansion should pull in more than just the zinc atom itself"
        );
    }

    #[test]
    fn whole_structure_selects_every_real_atom() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");

        let selection = Selection::whole_structure(&pdb);
        for hierarchy in pdb.atoms_with_hierarchy() {
            assert!(selection.contains(hierarchy.atom()));
        }
    }
}
