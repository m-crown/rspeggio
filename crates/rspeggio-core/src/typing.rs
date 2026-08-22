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

// The carbon at the center of a guanidinium/amidinium-style cation: bonded
// to exactly three nitrogens (arginine's CZ). The positive charge is
// delocalized across the whole C+3N group -- not localized on one atom --
// so this identifies the group by its shape, not by counting bonds on any
// single nitrogen alone.
fn is_guanidinium_carbon(atom: &CcdAtom, component: &CcdComponent) -> bool {
    if atom.element() != "C" {
        return false;
    }
    let neighbors = bonded_neighbors(atom, component);
    neighbors.len() == 3 && neighbors.iter().all(|(n, _)| n.element() == "N")
}

// True for a nitrogen that's either a fully-saturated (4-bond) ammonium,
// or one of the three nitrogens on a guanidinium carbon. Both are
// positively charged with no lone pair left to accept a hydrogen bond.
fn is_pos_ionisable_nitrogen(atom: &CcdAtom, component: &CcdComponent) -> bool {
    if atom.element() != "N" {
        return false;
    }
    let neighbors = bonded_neighbors(atom, component);
    neighbors.len() == 4
        || neighbors
            .iter()
            .any(|(n, _)| is_guanidinium_carbon(n, component))
}

// A carboxyl carbon: bonded to exactly two oxygens, one via a double bond
// and one via a single bond, where that single-bonded oxygen is terminal
// (its only other neighbor, if any, is a hydrogen -- not another carbon,
// which would make this an ester rather than a free/ionisable acid).
// Matches both the protonated (-COOH) and deprotonated (-COO-) forms,
// since "ionisable" means "capable of carrying charge", not "does here".
fn is_carboxyl_carbon(atom: &CcdAtom, component: &CcdComponent) -> bool {
    if atom.element() != "C" {
        return false;
    }
    let neighbors = bonded_neighbors(atom, component);
    let oxygens: Vec<_> = neighbors
        .iter()
        .filter(|(n, _)| n.element() == "O")
        .collect();
    if oxygens.len() != 2 {
        return false;
    }
    let has_double = oxygens
        .iter()
        .any(|(_, order)| **order == BondOrder::Double);
    let single_oxygen_is_terminal = oxygens
        .iter()
        .find(|(_, order)| **order != BondOrder::Double)
        .map(|(o, _)| {
            bonded_neighbors(o, component)
                .iter()
                .all(|(n, _)| n.element() == "H" || n.atom_id() == atom.atom_id())
        })
        .unwrap_or(false);
    has_double && single_oxygen_is_terminal
}

fn is_neg_ionisable_oxygen(atom: &CcdAtom, component: &CcdComponent) -> bool {
    if atom.element() != "O" {
        return false;
    }
    bonded_neighbors(atom, component)
        .iter()
        .any(|(n, _)| is_carboxyl_carbon(n, component))
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

    // Ionisable groups, identified structurally rather than by name --
    // guanidinium/ammonium cations (positive), free/deprotonatable
    // carboxyl groups (negative). These override the plain valence-based
    // acceptor rule below: a nitrogen with a technically-free lone pair by
    // bond count alone still isn't a real acceptor if that lone pair is
    // tied up carrying a formal positive charge.
    let pos_ionisable = is_pos_ionisable_nitrogen(atom, component)
        || (atom.element() == "C" && is_guanidinium_carbon(atom, component));
    let neg_ionisable = is_neg_ionisable_oxygen(atom, component);
    if pos_ionisable {
        bits |= AtomTypeBits::POS_IONISABLE;
    }
    if neg_ionisable {
        bits |= AtomTypeBits::NEG_IONISABLE;
    }

    // Donor/acceptor capacity from the free-form CCD chemistry graph, not
    // from whether any specific structure resolved a hydrogen -- most
    // X-ray structures have none at all. See the `donor_hydrogen_count_vs_terminus`
    // note: this is deliberately a boolean, not a hydrogen count, so
    // terminus status (which changes the count but not this boolean)
    // doesn't need to be known here.
    let bond_count = neighbors.len();
    let has_h_neighbor = neighbors.iter().any(|(other, _)| other.element() == "H");

    // Every oxygen has a lone pair free to accept, regardless of whether
    // it's also part of a negatively-ionisable carboxylate -- that charge
    // state doesn't consume the lone pair the way a cation's does.
    // Nitrogen has one too, unless all 4 bonding sites are used, or unless
    // it's positively ionisable (guanidinium's resonance ties up what
    // would otherwise look like a free lone pair by bond count alone).
    if atom.element() == "O" || (atom.element() == "N" && bond_count <= 3 && !pos_ionisable) {
        bits |= AtomTypeBits::HBOND_ACCEPTOR;
    }
    if matches!(atom.element(), "N" | "O" | "S") && has_h_neighbor {
        bits |= AtomTypeBits::HBOND_DONOR;
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
    fn alanines_backbone_amine_nitrogen_is_both_donor_and_acceptor() {
        let components = common_components();
        let ala = components.get("ALA").expect("ALA should be bundled");

        // N is bonded to CA, H, H2 -- 3 bonds, has an H neighbor.
        let n = ala
            .atoms()
            .iter()
            .find(|a| a.atom_id() == "N")
            .expect("ALA has a backbone N");

        let bits = type_atom(n, ala);
        assert!(bits.contains(AtomTypeBits::HBOND_DONOR));
        assert!(bits.contains(AtomTypeBits::HBOND_ACCEPTOR));
    }

    #[test]
    fn alanines_carbonyl_oxygen_is_acceptor_only() {
        let components = common_components();
        let ala = components.get("ALA").expect("ALA should be bundled");

        // O is double-bonded only to C -- no attached hydrogen.
        let o = ala
            .atoms()
            .iter()
            .find(|a| a.atom_id() == "O")
            .expect("ALA has a backbone carbonyl O");

        let bits = type_atom(o, ala);
        assert!(bits.contains(AtomTypeBits::HBOND_ACCEPTOR));
        assert!(!bits.contains(AtomTypeBits::HBOND_DONOR));
    }

    #[test]
    fn alanines_free_acid_hydroxyl_is_both_donor_and_acceptor() {
        let components = common_components();
        let ala = components.get("ALA").expect("ALA should be bundled");

        // OXT is bonded to C and HXT -- 2 bonds, has an H neighbor.
        let oxt = ala
            .atoms()
            .iter()
            .find(|a| a.atom_id() == "OXT")
            .expect("ALA has an OXT in its free-acid CCD form");

        let bits = type_atom(oxt, ala);
        assert!(bits.contains(AtomTypeBits::HBOND_DONOR));
        assert!(bits.contains(AtomTypeBits::HBOND_ACCEPTOR));
    }

    #[test]
    fn lysines_protonated_terminal_amine_is_donor_but_not_acceptor() {
        let components = common_components();
        let lys = components.get("LYS").expect("LYS should be bundled");

        // NZ is bonded to CE, HZ1, HZ2, HZ3 -- 4 bonds, no lone pair left
        // (this free-form CCD entry models the protonated NH3+ form).
        let nz = lys
            .atoms()
            .iter()
            .find(|a| a.atom_id() == "NZ")
            .expect("LYS has an NZ");

        let bits = type_atom(nz, lys);
        assert!(bits.contains(AtomTypeBits::HBOND_DONOR));
        assert!(
            !bits.contains(AtomTypeBits::HBOND_ACCEPTOR),
            "4 bonds leaves no lone pair to accept with"
        );
    }

    #[test]
    fn arginines_guanidinium_group_is_pos_ionisable_and_donor_but_never_acceptor() {
        let components = common_components();
        let arg = components.get("ARG").expect("ARG should be bundled");

        for id in ["CZ", "NE", "NH1", "NH2"] {
            let atom = arg
                .atoms()
                .iter()
                .find(|a| a.atom_id() == id)
                .unwrap_or_else(|| panic!("ARG should have a {id}"));
            let bits = type_atom(atom, arg);
            assert!(
                bits.contains(AtomTypeBits::POS_IONISABLE),
                "{id} should be part of the guanidinium group"
            );
            assert!(
                !bits.contains(AtomTypeBits::HBOND_ACCEPTOR),
                "{id}'s delocalized positive charge leaves no lone pair to accept with"
            );
        }
        // CZ itself has no attached hydrogen -- NE/NH1/NH2 do.
        let cz = arg.atoms().iter().find(|a| a.atom_id() == "CZ").unwrap();
        assert!(!type_atom(cz, arg).contains(AtomTypeBits::HBOND_DONOR));
        for id in ["NE", "NH1", "NH2"] {
            let atom = arg.atoms().iter().find(|a| a.atom_id() == id).unwrap();
            assert!(type_atom(atom, arg).contains(AtomTypeBits::HBOND_DONOR));
        }
    }

    #[test]
    fn lysines_nz_is_also_pos_ionisable() {
        let components = common_components();
        let lys = components.get("LYS").expect("LYS should be bundled");
        let nz = lys.atoms().iter().find(|a| a.atom_id() == "NZ").unwrap();

        assert!(type_atom(nz, lys).contains(AtomTypeBits::POS_IONISABLE));
    }

    #[test]
    fn aspartates_side_chain_carboxylate_is_neg_ionisable_and_still_an_acceptor() {
        let components = common_components();
        let asp = components.get("ASP").expect("ASP should be bundled");

        let od1 = asp.atoms().iter().find(|a| a.atom_id() == "OD1").unwrap();
        let od2 = asp.atoms().iter().find(|a| a.atom_id() == "OD2").unwrap();

        let od1_bits = type_atom(od1, asp);
        let od2_bits = type_atom(od2, asp);

        assert!(od1_bits.contains(AtomTypeBits::NEG_IONISABLE));
        assert!(od2_bits.contains(AtomTypeBits::NEG_IONISABLE));
        // negative ionisability doesn't consume the lone pair the way a
        // cation's positive charge does -- both stay acceptors, matching
        // arpeggio's own convention (ASPOD1/OD2 are both hbond-acceptor
        // and neg-ionisable simultaneously).
        assert!(od1_bits.contains(AtomTypeBits::HBOND_ACCEPTOR));
        assert!(od2_bits.contains(AtomTypeBits::HBOND_ACCEPTOR));
        // only OD2 carries the free-acid hydrogen in this CCD form.
        assert!(!od1_bits.contains(AtomTypeBits::HBOND_DONOR));
        assert!(od2_bits.contains(AtomTypeBits::HBOND_DONOR));
    }

    #[test]
    fn alanines_free_acid_backbone_is_also_neg_ionisable_unlike_arpeggios_own_table() {
        // Deliberate, known divergence from real arpeggio: its hardcoded
        // PROT_ATOM_TYPES only lists ASP/GLU side chains as neg-ionisable,
        // not any backbone terminus, even though a free residue's backbone
        // -COOH (this CCD entry's O/OXT) is structurally the same acid
        // group. Our rule identifies carboxyl groups generically by shape
        // rather than by a curated per-residue list, so it correctly (if
        // more broadly than the oracle) flags this one too.
        let components = common_components();
        let ala = components.get("ALA").expect("ALA should be bundled");
        let o = ala.atoms().iter().find(|a| a.atom_id() == "O").unwrap();
        let oxt = ala.atoms().iter().find(|a| a.atom_id() == "OXT").unwrap();

        assert!(type_atom(o, ala).contains(AtomTypeBits::NEG_IONISABLE));
        assert!(type_atom(oxt, ala).contains(AtomTypeBits::NEG_IONISABLE));
    }

    #[test]
    fn a_plain_carbon_is_neither_donor_nor_acceptor() {
        let components = common_components();
        let ala = components.get("ALA").expect("ALA should be bundled");

        let cb = ala
            .atoms()
            .iter()
            .find(|a| a.atom_id() == "CB")
            .expect("ALA has a CB");

        let bits = type_atom(cb, ala);
        assert!(!bits.contains(AtomTypeBits::HBOND_DONOR));
        assert!(!bits.contains(AtomTypeBits::HBOND_ACCEPTOR));
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
