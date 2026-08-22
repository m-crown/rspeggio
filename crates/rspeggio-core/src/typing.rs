// crates/rspeggio-core/src/typing.rs
//
// Hand-rolled atom typing (M4): classifies a CCD atom's chemical role
// (hbond donor/acceptor, ionisable, hydrophobe, carbonyl, aromatic, ...)
// directly from its element, bonds, and valence -- no SMARTS matcher, no
// chemistry-library FFI, per decisions 01/02.
//
// Typing is a property of the *component's chemistry graph*, not of any
// particular structure instance: the same MET residue types identically
// everywhere it appears. So this operates on a `CcdAtom` + its owning
// `CcdComponent` (for bond lookups), not on a `JoinedAtom` -- the join's
// job was only to find *which* CCD atom corresponds to a structure atom;
// typing itself needs nothing structure-specific.

use bitflags::bitflags;
use rspeggio_ccd::component::{BondOrder, CcdAtom, CcdComponent};

bitflags! {
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub struct AtomTypeBits: u16 {
        const HBOND_DONOR         = 1 << 0;
        const HBOND_ACCEPTOR      = 1 << 1;
        const WEAK_HBOND_DONOR    = 1 << 2;
        const WEAK_HBOND_ACCEPTOR = 1 << 3;
        const POS_IONISABLE       = 1 << 4;
        const NEG_IONISABLE       = 1 << 5;
        const HYDROPHOBE          = 1 << 6;
        const AROMATIC            = 1 << 7;
        const CARBONYL_OXYGEN     = 1 << 8;
        const CARBONYL_CARBON     = 1 << 9;
        const METAL               = 1 << 10;
    }
}

// Elements treated as metals for typing purposes. Covers what's in the
// bundled common set (Na, Mg, Ca, K) plus the other biologically common
// ones the fixture corpus references (Zn, Fe) -- not an exhaustive
// periodic-table list, extend as real structures need more.
const METAL_ELEMENTS: &[&str] = &["NA", "MG", "CA", "K", "ZN", "FE", "MN", "CU", "CO", "NI"];

// Every other atom this atom is bonded to within its component, paired
// with the bond order connecting them. Looked up by name since that's how
// `CcdBond` records its endpoints -- there's no atom-id-indexed adjacency
// list yet, so this is a linear scan; fine at CCD-component scale (single
// digits to low hundreds of atoms), worth revisiting only if typing ever
// needs to run on much larger components.
fn bonded_neighbors<'a>(
    atom: &CcdAtom,
    component: &'a CcdComponent,
) -> Vec<(&'a CcdAtom, &'a BondOrder)> {
    component
        .bonds()
        .iter()
        .filter_map(|bond| {
            let other_id = if bond.atom_id_1() == atom.atom_id() {
                Some(bond.atom_id_2())
            } else if bond.atom_id_2() == atom.atom_id() {
                Some(bond.atom_id_1())
            } else {
                None
            }?;
            component
                .atoms()
                .iter()
                .find(|a| a.atom_id() == other_id)
                .map(|other| (other, bond.order()))
        })
        .collect()
}

pub fn type_atom(atom: &CcdAtom, component: &CcdComponent) -> AtomTypeBits {
    let mut bits = AtomTypeBits::empty();

    if atom.aromatic() {
        bits |= AtomTypeBits::AROMATIC;
    }

    if METAL_ELEMENTS.contains(&atom.element().to_uppercase().as_str()) {
        bits |= AtomTypeBits::METAL;
    }

    let neighbors = bonded_neighbors(atom, component);

    // A carbonyl is the C=O pair itself: this atom is the O double-bonded
    // to a C, or the C double-bonded to an O -- not "any oxygen near a
    // carbon", which would also wrongly catch e.g. a C-OH hydroxyl.
    if atom.element() == "O"
        && neighbors
            .iter()
            .any(|(other, order)| other.element() == "C" && **order == BondOrder::Double)
    {
        bits |= AtomTypeBits::CARBONYL_OXYGEN;
    }
    if atom.element() == "C"
        && neighbors
            .iter()
            .any(|(other, order)| other.element() == "O" && **order == BondOrder::Double)
    {
        bits |= AtomTypeBits::CARBONYL_CARBON;
    }

    // A hydrophobe is a carbon with no heteroatom neighbors at all -- every
    // bond goes to another carbon or to hydrogen. This is a property of
    // the whole neighbor set, not any single bond, so (unlike carbonyl)
    // it's `all`, not `any`. A carbon with zero bonds (shouldn't occur in
    // real CCD data, but `all` on an empty iterator is vacuously true)
    // wouldn't count -- guarded explicitly rather than relying on that.
    if atom.element() == "C"
        && !neighbors.is_empty()
        && neighbors
            .iter()
            .all(|(other, _)| matches!(other.element(), "C" | "H"))
    {
        bits |= AtomTypeBits::HYDROPHOBE;
    }

    bits
}

#[cfg(test)]
mod tests {
    use super::*;
    use rspeggio_ccd::common::common_components;

    #[test]
    fn phenylalanines_ring_carbons_are_typed_aromatic() {
        let components = common_components();
        let phe = components.get("PHE").expect("PHE should be bundled");

        let cz = phe
            .atoms()
            .iter()
            .find(|a| a.atom_id() == "CZ")
            .expect("PHE has a CZ ring atom");

        let bits = type_atom(cz, phe);

        assert!(bits.contains(AtomTypeBits::AROMATIC));
        assert!(!bits.contains(AtomTypeBits::METAL));
    }

    #[test]
    fn phenylalanines_backbone_nitrogen_is_not_aromatic() {
        let components = common_components();
        let phe = components.get("PHE").expect("PHE should be bundled");

        let n = phe
            .atoms()
            .iter()
            .find(|a| a.atom_id() == "N")
            .expect("PHE has a backbone N");

        let bits = type_atom(n, phe);

        assert!(!bits.contains(AtomTypeBits::AROMATIC));
    }

    #[test]
    fn bundled_calcium_ion_is_typed_metal() {
        let components = common_components();
        let ca_ion = components.get("CA").expect("CA should be bundled");

        let ca_atom = ca_ion
            .atoms()
            .iter()
            .find(|a| a.atom_id() == "CA")
            .expect("the CA component has one atom named CA");

        let bits = type_atom(ca_atom, ca_ion);

        assert!(bits.contains(AtomTypeBits::METAL));
        assert!(!bits.contains(AtomTypeBits::AROMATIC));
    }

    #[test]
    fn alanines_backbone_carbonyl_is_typed_on_both_atoms() {
        let components = common_components();
        let ala = components.get("ALA").expect("ALA should be bundled");

        let c = ala
            .atoms()
            .iter()
            .find(|a| a.atom_id() == "C")
            .expect("ALA has a backbone C");
        let o = ala
            .atoms()
            .iter()
            .find(|a| a.atom_id() == "O")
            .expect("ALA has a backbone carbonyl O");

        assert!(type_atom(c, ala).contains(AtomTypeBits::CARBONYL_CARBON));
        assert!(type_atom(o, ala).contains(AtomTypeBits::CARBONYL_OXYGEN));
    }

    #[test]
    fn alanines_hydroxyl_oxygen_is_not_a_carbonyl() {
        let components = common_components();
        let ala = components.get("ALA").expect("ALA should be bundled");

        // OXT is singly bonded to C (the free acid's -OH), not a carbonyl,
        // even though it's bonded to the same carbon that has one.
        let oxt = ala
            .atoms()
            .iter()
            .find(|a| a.atom_id() == "OXT")
            .expect("ALA has an OXT in its free-acid CCD form");

        let bits = type_atom(oxt, ala);

        assert!(!bits.contains(AtomTypeBits::CARBONYL_OXYGEN));
        assert!(!bits.contains(AtomTypeBits::CARBONYL_CARBON));
    }

    #[test]
    fn alanines_alpha_carbon_is_not_a_carbonyl_carbon() {
        let components = common_components();
        let ala = components.get("ALA").expect("ALA should be bundled");

        // CA is bonded to N, C, CB, HA -- none of them an O, let alone a
        // double bond to one.
        let ca = ala
            .atoms()
            .iter()
            .find(|a| a.atom_id() == "CA")
            .expect("ALA has a CA");

        assert!(!type_atom(ca, ala).contains(AtomTypeBits::CARBONYL_CARBON));
    }

    #[test]
    fn alanines_methyl_carbon_is_hydrophobe() {
        let components = common_components();
        let ala = components.get("ALA").expect("ALA should be bundled");

        // CB is bonded only to CA, HB1, HB2, HB3 -- carbon and hydrogens.
        let cb = ala
            .atoms()
            .iter()
            .find(|a| a.atom_id() == "CB")
            .expect("ALA has a CB");

        assert!(type_atom(cb, ala).contains(AtomTypeBits::HYDROPHOBE));
    }

    #[test]
    fn alanines_alpha_carbon_is_not_hydrophobe_because_it_touches_nitrogen() {
        let components = common_components();
        let ala = components.get("ALA").expect("ALA should be bundled");

        // CA is bonded to N, C, CB, HA -- N is a heteroatom neighbor.
        let ca = ala
            .atoms()
            .iter()
            .find(|a| a.atom_id() == "CA")
            .expect("ALA has a CA");

        assert!(!type_atom(ca, ala).contains(AtomTypeBits::HYDROPHOBE));
    }

    #[test]
    fn alanines_carbonyl_carbon_is_not_hydrophobe_because_it_touches_oxygen() {
        let components = common_components();
        let ala = components.get("ALA").expect("ALA should be bundled");

        let c = ala
            .atoms()
            .iter()
            .find(|a| a.atom_id() == "C")
            .expect("ALA has a backbone C");

        assert!(!type_atom(c, ala).contains(AtomTypeBits::HYDROPHOBE));
    }

    #[test]
    fn phenylalanines_ring_carbons_are_also_hydrophobe() {
        let components = common_components();
        let phe = components.get("PHE").expect("PHE should be bundled");

        // CZ (ring apex) is bonded only to CE1, CE2, HZ -- an aromatic
        // carbon can be both AROMATIC and HYDROPHOBE, they're independent
        // bits, not mutually exclusive categories.
        let cz = phe
            .atoms()
            .iter()
            .find(|a| a.atom_id() == "CZ")
            .expect("PHE has a CZ");

        let bits = type_atom(cz, phe);
        assert!(bits.contains(AtomTypeBits::AROMATIC));
        assert!(bits.contains(AtomTypeBits::HYDROPHOBE));
    }

    #[test]
    fn a_plain_carbon_is_neither_aromatic_nor_metal() {
        let components = common_components();
        let ala = components.get("ALA").expect("ALA should be bundled");

        let cb = ala
            .atoms()
            .iter()
            .find(|a| a.atom_id() == "CB")
            .expect("ALA has a CB atom");

        let bits = type_atom(cb, ala);

        // CB is genuinely HYDROPHOBE (see alanines_methyl_carbon_is_hydrophobe)
        // -- this test only cares that the unrelated categories stay clear.
        assert!(!bits.contains(AtomTypeBits::AROMATIC));
        assert!(!bits.contains(AtomTypeBits::METAL));
        assert!(!bits.contains(AtomTypeBits::CARBONYL_CARBON));
        assert!(!bits.contains(AtomTypeBits::CARBONYL_OXYGEN));
    }
}
