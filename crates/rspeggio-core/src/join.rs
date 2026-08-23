// crates/rspeggio-core/src/join.rs
//
// Cross-references pdbtbx's parsed structure atoms (pure geometry) against
// real chemistry from CCD `CcdComponent`s (bonds, aromaticity, from
// rspeggio-ccd). The join key is `(residue comp_id, atom name)`: a
// residue's name selects the component, and each atom is matched to the
// CCD atom of the same name within it.
//
// This operates directly on pdbtbx's own hierarchy types
// (`AtomConformerResidueChainModel`) rather than a separate flattened
// `Structure` abstraction -- pdbtbx's own model already carries everything
// needed (atom/residue/chain context, plus rstar-based spatial queries via
// its `rstar` feature), so a parallel representation would just be
// redundant upkeep.
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

use pdbtbx::AtomConformerResidueChainModel;
use pdbtbx::{ContainsAtomConformer, ContainsAtomConformerResidue, PDB};
use rspeggio_ccd::component::{CcdAtom, CcdComponent};
use std::collections::HashMap;

#[derive(Clone)]
pub struct JoinedAtom<'a> {
    pub hierarchy: AtomConformerResidueChainModel<'a>,
    pub component: Option<&'a CcdComponent>,
    pub ccd_atom: Option<&'a CcdAtom>,
}

pub fn join_atom<'a>(
    hierarchy: AtomConformerResidueChainModel<'a>,
    components: &'a HashMap<String, CcdComponent>,
) -> JoinedAtom<'a> {
    let comp_id = hierarchy.residue().name().unwrap_or_default();
    let component = components.get(comp_id);
    let ccd_atom = component.and_then(|c| {
        c.atoms()
            .iter()
            .find(|a| a.atom_id() == hierarchy.atom().name())
    });
    JoinedAtom {
        hierarchy,
        component,
        ccd_atom,
    }
}

// Joins every atom in a structure. Order matches `pdb.atoms_with_hierarchy()`.
pub fn join_all<'a>(
    pdb: &'a PDB,
    components: &'a HashMap<String, CcdComponent>,
) -> Vec<JoinedAtom<'a>> {
    pdb.atoms_with_hierarchy()
        .map(|hierarchy| join_atom(hierarchy, components))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rspeggio_ccd::common::common_components;

    #[test]
    fn joins_a_known_residues_atoms_to_their_ccd_chemistry() {
        let (pdb, _errors) = pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif")
            .expect("1UBQ fixture should parse");
        let components = common_components();

        let joined = join_all(&pdb, &components);

        let ca = joined
            .iter()
            .find(|j| j.hierarchy.residue().id().0 == 1 && j.hierarchy.atom().name() == "CA")
            .expect("MET1 CA should be in the joined atom list");

        assert!(ca.component.is_some(), "MET is in the common bundle");
        let ccd_ca = ca.ccd_atom.expect("CA should match a CCD atom by name");
        assert_eq!(ccd_ca.element(), "C");
        assert!(!ccd_ca.aromatic());
    }

    #[test]
    fn joins_every_atom_in_a_real_structure() {
        let (pdb, _errors) = pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif")
            .expect("1UBQ fixture should parse");
        let components = common_components();

        let joined = join_all(&pdb, &components);

        assert_eq!(joined.len(), pdb.atom_count());
        assert!(
            joined.iter().all(|j| j.component.is_some()),
            "1UBQ is amino acids + water only, all in the common bundle"
        );
        assert!(
            joined.iter().all(|j| j.ccd_atom.is_some()),
            "every real atom in this protein-only fixture should match its CCD entry by name"
        );
    }

    #[test]
    fn a_residue_with_no_matching_component_joins_with_none_rather_than_failing() {
        let (pdb, _errors) = pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif")
            .expect("1UBQ fixture should parse");
        // Deliberately empty -- simulates every residue being unknown.
        let components: HashMap<String, CcdComponent> = HashMap::new();

        let joined = join_all(&pdb, &components);

        assert!(!joined.is_empty());
        assert!(joined.iter().all(|j| j.component.is_none()));
        assert!(joined.iter().all(|j| j.ccd_atom.is_none()));
    }

    #[test]
    fn an_atom_absent_from_an_otherwise_known_component_joins_with_none() {
        let (pdb, _errors) = pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif")
            .expect("1UBQ fixture should parse");

        // 1UBQ has no hydrogens at all (a normal X-ray structure), so
        // every real atom it contains matches its CCD entry by name --
        // there's no naturally-occurring "atom missing from a known
        // component" case in this fixture. Build one honestly instead:
        // take MET's real bundled atoms and bonds, but drop CA from the
        // atom list, so the component is still real chemistry data, just
        // deliberately incomplete for this one atom.
        let mut components = common_components();
        let met = components.remove("MET").expect("MET should be bundled");
        let atoms: Vec<CcdAtom> = met
            .atoms()
            .iter()
            .filter(|a| a.atom_id() != "CA")
            .cloned()
            .collect();
        let bonds: Vec<_> = met.bonds().to_vec();
        components.insert("MET".to_string(), CcdComponent::new(atoms, bonds));

        let joined = join_all(&pdb, &components);

        let ca = joined
            .iter()
            .find(|j| j.hierarchy.residue().id().0 == 1 && j.hierarchy.atom().name() == "CA")
            .expect("MET1 CA should be present in the structure");
        assert!(ca.component.is_some(), "MET itself is still known");
        assert!(
            ca.ccd_atom.is_none(),
            "CA was deliberately removed from this component's atom list"
        );

        let n = joined
            .iter()
            .find(|j| j.hierarchy.residue().id().0 == 1 && j.hierarchy.atom().name() == "N")
            .expect("MET1 N should be present");
        assert!(
            n.ccd_atom.is_some(),
            "other atoms in the same residue still match"
        );
    }
}
