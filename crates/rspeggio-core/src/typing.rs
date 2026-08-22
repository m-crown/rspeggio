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
use rspeggio_ccd::component::{CcdAtom, CcdComponent};

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

pub fn type_atom(atom: &CcdAtom, _component: &CcdComponent) -> AtomTypeBits {
    let mut bits = AtomTypeBits::empty();

    if atom.aromatic() {
        bits |= AtomTypeBits::AROMATIC;
    }

    if METAL_ELEMENTS.contains(&atom.element().to_uppercase().as_str()) {
        bits |= AtomTypeBits::METAL;
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
    fn a_plain_carbon_is_neither_aromatic_nor_metal() {
        let components = common_components();
        let ala = components.get("ALA").expect("ALA should be bundled");

        let cb = ala
            .atoms()
            .iter()
            .find(|a| a.atom_id() == "CB")
            .expect("ALA has a CB atom");

        let bits = type_atom(cb, ala);

        assert!(bits.is_empty());
    }
}
