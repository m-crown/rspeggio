// crates/rspeggio-core/src/join.rs
//
// Cross-references a `Structure`'s residues/atoms (pure geometry, from
// pdbtbx) against real chemistry from CCD `CcdComponent`s (bonds,
// aromaticity, from rspeggio-ccd). The join key is `(comp_id, atom name)`:
// a structure residue's `comp_id` selects the component, and each atom is
// matched to the CCD atom of the same name within it.
//
// Two independent things can be missing, and both are represented as
// `None` rather than failing the whole join -- a structure legitimately
// mixes well-known residues with occasional gaps:
//   - the residue's `comp_id` has no matching component at all (e.g. an
//     uncommon ligand not in the bundled set -- the disk-cache/live-fetch
//     tier that will fill this gap is separate, later work)
//   - the component exists, but this particular atom isn't in it (e.g. a
//     hydrogen present in the structure but absent from CCD ideal
//     coordinates, which omit hydrogens for many entries)

use rspeggio_ccd::component::{CcdAtom, CcdComponent};
use std::collections::HashMap;

use crate::structure::{Structure, StructureAtom, StructureResidue};

pub struct JoinedAtom<'a> {
    pub structure_atom: &'a StructureAtom,
    pub ccd_atom: Option<&'a CcdAtom>,
}

pub struct JoinedResidue<'a> {
    pub residue: &'a StructureResidue,
    pub component: Option<&'a CcdComponent>,
    pub atoms: Vec<JoinedAtom<'a>>,
}

pub fn join_residue<'a>(
    residue: &'a StructureResidue,
    components: &'a HashMap<String, CcdComponent>,
) -> JoinedResidue<'a> {
    let component = components.get(&residue.comp_id);

    let atoms = residue
        .atoms
        .iter()
        .map(|structure_atom| {
            let ccd_atom = component.and_then(|c| {
                c.atoms()
                    .iter()
                    .find(|a| a.atom_id() == structure_atom.name)
            });
            JoinedAtom {
                structure_atom,
                ccd_atom,
            }
        })
        .collect();

    JoinedResidue {
        residue,
        component,
        atoms,
    }
}

// Joins every residue in a structure. Order matches `structure.residues`.
pub fn join_structure<'a>(
    structure: &'a Structure,
    components: &'a HashMap<String, CcdComponent>,
) -> Vec<JoinedResidue<'a>> {
    structure
        .residues
        .iter()
        .map(|residue| join_residue(residue, components))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::load_structure;
    use rspeggio_ccd::common::common_components;

    #[test]
    fn joins_a_known_residues_atoms_to_their_ccd_chemistry() {
        let structure =
            load_structure("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");
        let components = common_components();

        let met1 = structure
            .residues
            .iter()
            .find(|r| r.comp_id == "MET" && r.seq_id == 1)
            .expect("MET1 should be present");

        let joined = join_residue(met1, &components);

        assert!(joined.component.is_some(), "MET is in the common bundle");

        let ca = joined
            .atoms
            .iter()
            .find(|a| a.structure_atom.name == "CA")
            .expect("CA should be in the joined atom list");
        let ccd_ca = ca.ccd_atom.expect("CA should match a CCD atom by name");
        assert_eq!(ccd_ca.element(), "C");
        assert!(!ccd_ca.aromatic());
    }

    #[test]
    fn joins_every_residue_in_a_real_structure() {
        let structure =
            load_structure("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");
        let components = common_components();

        let joined = join_structure(&structure, &components);

        assert_eq!(joined.len(), structure.residues.len());
        assert!(
            joined.iter().all(|r| r.component.is_some()),
            "1UBQ is amino acids + water only, all in the common bundle"
        );

        let total_atoms: usize = joined.iter().map(|r| r.atoms.len()).sum();
        let matched_atoms = joined
            .iter()
            .flat_map(|r| &r.atoms)
            .filter(|a| a.ccd_atom.is_some())
            .count();
        assert_eq!(
            total_atoms, matched_atoms,
            "every real atom in this protein-only fixture should match its CCD entry by name"
        );
    }

    #[test]
    fn a_residue_with_no_matching_component_joins_with_none_rather_than_failing() {
        let residue = StructureResidue {
            comp_id: "ZZZ9".to_string(),
            chain_id: "A".to_string(),
            seq_id: 1,
            atoms: vec![StructureAtom {
                serial: 1,
                name: "X1".to_string(),
                element: "X".to_string(),
                pos: (0.0, 0.0, 0.0),
            }],
        };
        let components = common_components();

        let joined = join_residue(&residue, &components);

        assert!(joined.component.is_none());
        assert_eq!(joined.atoms.len(), 1);
        assert!(joined.atoms[0].ccd_atom.is_none());
    }

    #[test]
    fn a_structure_atom_absent_from_an_otherwise_known_component_joins_with_none() {
        let residue = StructureResidue {
            comp_id: "MET".to_string(),
            chain_id: "A".to_string(),
            seq_id: 1,
            atoms: vec![StructureAtom {
                serial: 1,
                name: "NOT_A_REAL_ATOM_NAME".to_string(),
                element: "C".to_string(),
                pos: (0.0, 0.0, 0.0),
            }],
        };
        let components = common_components();

        let joined = join_residue(&residue, &components);

        assert!(joined.component.is_some(), "MET itself is known");
        assert!(
            joined.atoms[0].ccd_atom.is_none(),
            "this specific atom name isn't in MET's CCD entry"
        );
    }
}
