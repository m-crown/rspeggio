// crates/rspeggio-core/src/features.rs
//
// The feature half of M5's atom-atom contact classification (SIFt
// positions 5-13: HBOND, WEAK_HBOND, XBOND, IONIC, METAL_COMPLEX,
// AROMATIC, HYDROPHOBIC, CARBONYL), built on top of `contacts.rs`'s
// distance-category half. This is where M4's atom typing (`AtomTypeBits`)
// finally gets consumed: each contact's two atoms are joined to their CCD
// chemistry, typed, then checked against the distance thresholds already
// ported in `config.rs`.
//
// IONIC, HYDROPHOBIC, CARBONYL, METAL_COMPLEX, AROMATIC need only distance
// + typing (AROMATIC here is real pdbe-arpeggio's simple atom-atom version,
// `interactions.py:916` -- both atoms typed aromatic within 4.0A; the
// separate ring-plane/centroid geometry that M6 also needed lives in
// `rings.rs` as its own plane-plane contact type, not an atom-pair
// feature). HBOND/WEAK_HBOND need a real donor-hydrogen position to check
// the donor-H...acceptor angle -- this module uses `hydrogenate` to place
// one wherever the donor's CCD shape is covered (real or resolved atoms
// take priority over a placed one), and falls back to the plan's own
// angle-free "polar distance" convention (config::HBOND/WEAK_HBOND's
// `polar_distance`/`weak_polar_distance`) wherever it isn't -- most
// notably backbone amide N-H, which needs the previous residue's real C
// (a cross-residue lookup this module doesn't attempt), and water, which
// `hydrogenate`'s own module doc explains needs a different mechanism
// entirely. XBOND (halogen bond donor/acceptor) isn't typed yet, so it's
// still not covered here.

use crate::config::{self, FeatureBits};
use crate::contacts::{euclidean_distance, Contact};
use crate::hydrogenate::{self, Hybridization};
use crate::join::join_atom;
use crate::typing::{self, AtomTypeBits};
use pdbtbx::{
    AtomConformerResidueChainModel, ContainsAtomConformer, ContainsAtomConformerResidue, Residue,
};
use rspeggio_ccd::component::{BondOrder, CcdAtom, CcdComponent};
use std::collections::HashMap;

type Point = (f64, f64, f64);

struct TypedAtom<'a> {
    hierarchy: AtomConformerResidueChainModel<'a>,
    ccd_atom: Option<&'a CcdAtom>,
    component: Option<&'a CcdComponent>,
    bits: AtomTypeBits,
}

fn typed_atom<'a>(
    hierarchy: AtomConformerResidueChainModel<'a>,
    components: &'a HashMap<String, CcdComponent>,
) -> TypedAtom<'a> {
    let joined = join_atom(hierarchy.clone(), components);
    let bits = match (joined.ccd_atom, joined.component) {
        (Some(a), Some(c)) => typing::type_atom(a, c),
        _ => AtomTypeBits::empty(),
    };
    TypedAtom {
        hierarchy,
        ccd_atom: joined.ccd_atom,
        component: joined.component,
        bits,
    }
}

// True if the two type-bit sets satisfy the role pair in either atom
// order -- a contact's two atoms aren't ordered by role (donor/acceptor,
// carbon/oxygen, ...), so every one of these rules needs to check both
// ways rather than assuming which atom is which.
fn either_order(
    bits_1: AtomTypeBits,
    bits_2: AtomTypeBits,
    role_1: AtomTypeBits,
    role_2: AtomTypeBits,
) -> bool {
    (bits_1.contains(role_1) && bits_2.contains(role_2))
        || (bits_2.contains(role_1) && bits_1.contains(role_2))
}

// The real position of some other real, resolved neighbor of `atom` in
// the CCD bond graph (excluding `exclude_id`) -- used to orient a planar
// amine's shared plane (see `hydrogenate::place_planar_amine_hydrogens`).
fn reference_neighbor_position(
    atom: &CcdAtom,
    exclude_id: &str,
    component: &CcdComponent,
    residue: &Residue,
) -> Option<Point> {
    typing::bonded_neighbors(atom, component)
        .into_iter()
        .filter(|(n, _)| n.element() != "H" && n.atom_id() != exclude_id)
        .find_map(|(n, _)| {
            residue
                .atoms()
                .find(|a| a.name() == n.atom_id())
                .map(|a| a.pos())
        })
}

// Finds the real 3D positions of every hydrogen a donor atom's CCD entry
// says it has: the real, resolved atom when the structure actually has
// one by that name, otherwise placed analytically via `hydrogenate` using
// only real neighbor positions already resolved in the *same residue*.
// Returns whatever it could resolve -- possibly fewer positions than the
// CCD lists, or none at all for shapes not covered (backbone N-H, water).
fn donor_hydrogen_positions(
    donor_ccd: &CcdAtom,
    donor_pos: Point,
    component: &CcdComponent,
    residue: &Residue,
) -> Vec<Point> {
    let neighbors = typing::bonded_neighbors(donor_ccd, component);
    let h_names: Vec<&str> = neighbors
        .iter()
        .filter(|(a, _)| a.element() == "H")
        .map(|(a, _)| a.atom_id())
        .collect();
    if h_names.is_empty() {
        return Vec::new();
    }

    let real_pos = |name: &str| -> Option<Point> {
        residue.atoms().find(|a| a.name() == name).map(|a| a.pos())
    };

    let mut positions = Vec::new();
    let mut missing = Vec::new();
    for name in &h_names {
        match real_pos(name) {
            Some(p) => positions.push(p),
            None => missing.push(*name),
        }
    }
    if missing.is_empty() {
        return positions;
    }

    let heavy_neighbors: Vec<&CcdAtom> = neighbors
        .iter()
        .filter(|(a, _)| a.element() != "H")
        .map(|(a, _)| *a)
        .collect();
    let heavy_positions: Vec<Point> = heavy_neighbors
        .iter()
        .filter_map(|n| real_pos(n.atom_id()))
        .collect();
    if heavy_positions.len() != heavy_neighbors.len() {
        // some heavy neighbor isn't resolved in this structure either --
        // can't reliably place anything more, return whatever real H's
        // were found (possibly none).
        return positions;
    }

    let bond_length = config::covalent_radius(donor_ccd.element()).unwrap_or(0.0)
        + config::covalent_radius("H").unwrap_or(0.0);
    let is_sp2 = donor_ccd.aromatic()
        || neighbors
            .iter()
            .any(|(_, order)| **order == BondOrder::Double);
    let hybridization = if is_sp2 {
        Hybridization::Sp2
    } else {
        Hybridization::Sp3
    };

    match (heavy_positions.len(), missing.len()) {
        (1, 1) => positions.push(hydrogenate::place_arbitrary_single_neighbor_substituent(
            donor_pos,
            heavy_positions[0],
            hybridization,
            bond_length,
        )),
        (1, 2) => {
            if let Some(reference) = reference_neighbor_position(
                heavy_neighbors[0],
                donor_ccd.atom_id(),
                component,
                residue,
            ) {
                let (h1, h2) = hydrogenate::place_planar_amine_hydrogens(
                    donor_pos,
                    heavy_positions[0],
                    reference,
                    bond_length,
                );
                positions.push(h1);
                positions.push(h2);
            }
        }
        (1, 3) => {
            let (h1, h2, h3) =
                hydrogenate::place_methyl_hydrogens(donor_pos, heavy_positions[0], bond_length);
            positions.push(h1);
            positions.push(h2);
            positions.push(h3);
        }
        (2, 1) => positions.push(hydrogenate::place_planar_two_neighbor_substituent(
            donor_pos,
            heavy_positions[0],
            heavy_positions[1],
            bond_length,
        )),
        (2, 2) => {
            let (h1, h2) = hydrogenate::place_tetrahedral_pair_from_two_neighbors(
                donor_pos,
                heavy_positions[0],
                heavy_positions[1],
                bond_length,
            );
            positions.push(h1);
            positions.push(h2);
        }
        (3, 1) => positions.push(hydrogenate::place_fourth_tetrahedral_substituent(
            donor_pos,
            heavy_positions[0],
            heavy_positions[1],
            heavy_positions[2],
            bond_length,
        )),
        _ => {} // unsupported shape (e.g. 0 known neighbors -- water)
    }

    positions
}

// `Some(true/false)` if real donor-H geometry could be checked (a real
// answer either way); `None` if no donor-H position could be placed at
// all, meaning the caller should fall back to the angle-free
// polar-distance convention instead.
fn donor_reaches_acceptor_by_angle(
    donor: &TypedAtom,
    acceptor: &TypedAtom,
    angle_threshold_deg: f64,
) -> Option<bool> {
    let donor_ccd = donor.ccd_atom?;
    let donor_component = donor.component?;
    let acceptor_ccd = acceptor.ccd_atom?;
    let donor_pos = donor.hierarchy.atom().pos();
    let acceptor_pos = acceptor.hierarchy.atom().pos();

    let h_positions = donor_hydrogen_positions(
        donor_ccd,
        donor_pos,
        donor_component,
        donor.hierarchy.residue(),
    );
    if h_positions.is_empty() {
        return None;
    }

    let max_h_dist = config::vdw_radius("H").unwrap_or(0.0)
        + config::vdw_radius(acceptor_ccd.element()).unwrap_or(0.0)
        + config::VDW_COMP_FACTOR;

    Some(h_positions.iter().any(|&h_pos| {
        euclidean_distance(h_pos, acceptor_pos) <= max_h_dist
            && hydrogenate::angle_degrees(donor_pos, h_pos, acceptor_pos) >= angle_threshold_deg
    }))
}

#[allow(clippy::too_many_arguments)]
fn classify_hbond_like(
    a1: &TypedAtom,
    a2: &TypedAtom,
    distance: f64,
    donor_bit: AtomTypeBits,
    acceptor_bit: AtomTypeBits,
    heavy_distance_threshold: f64,
    polar_distance_threshold: f64,
    angle_threshold_deg: f64,
) -> bool {
    if distance > heavy_distance_threshold {
        return false;
    }

    for (donor, acceptor) in [(a1, a2), (a2, a1)] {
        if !donor.bits.contains(donor_bit) || !acceptor.bits.contains(acceptor_bit) {
            continue;
        }
        match donor_reaches_acceptor_by_angle(donor, acceptor, angle_threshold_deg) {
            Some(true) => return true,
            Some(false) => {} // real geometry checked and failed this direction
            None => {
                // no donor-H could be placed -- fall back to the
                // angle-free polar-distance convention for this pair.
                if distance <= polar_distance_threshold {
                    return true;
                }
            }
        }
    }
    false
}

// Classifies every distance-and-typing feature contact for one atom pair
// (all of them except XBOND, which isn't typed yet).
pub fn classify_features(
    contact: &Contact,
    components: &HashMap<String, CcdComponent>,
) -> FeatureBits {
    let a1 = typed_atom(contact.atom_1.clone(), components);
    let a2 = typed_atom(contact.atom_2.clone(), components);

    let mut features = FeatureBits::empty();

    if contact.distance <= config::IONIC.distance
        && either_order(
            a1.bits,
            a2.bits,
            AtomTypeBits::POS_IONISABLE,
            AtomTypeBits::NEG_IONISABLE,
        )
    {
        features |= FeatureBits::IONIC;
    }

    if contact.distance <= config::HYDROPHOBIC.distance
        && a1.bits.contains(AtomTypeBits::HYDROPHOBE)
        && a2.bits.contains(AtomTypeBits::HYDROPHOBE)
    {
        features |= FeatureBits::HYDROPHOBIC;
    }

    if contact.distance <= config::AROMATIC.distance
        && a1.bits.contains(AtomTypeBits::AROMATIC)
        && a2.bits.contains(AtomTypeBits::AROMATIC)
    {
        features |= FeatureBits::AROMATIC;
    }

    if contact.distance <= config::CARBONYL.distance
        && either_order(
            a1.bits,
            a2.bits,
            AtomTypeBits::CARBONYL_CARBON,
            AtomTypeBits::CARBONYL_OXYGEN,
        )
    {
        features |= FeatureBits::CARBONYL;
    }

    if contact.distance <= config::METAL.distance
        && either_order(
            a1.bits,
            a2.bits,
            AtomTypeBits::METAL,
            AtomTypeBits::HBOND_ACCEPTOR,
        )
    {
        features |= FeatureBits::METAL_COMPLEX;
    }

    if classify_hbond_like(
        &a1,
        &a2,
        contact.distance,
        AtomTypeBits::HBOND_DONOR,
        AtomTypeBits::HBOND_ACCEPTOR,
        config::HBOND.distance,
        config::HBOND.polar_distance,
        config::HBOND.angle_degrees,
    ) {
        features |= FeatureBits::HBOND;
    }

    if classify_hbond_like(
        &a1,
        &a2,
        contact.distance,
        AtomTypeBits::WEAK_HBOND_DONOR,
        AtomTypeBits::WEAK_HBOND_ACCEPTOR,
        config::WEAK_HBOND.distance,
        config::WEAK_HBOND.weak_polar_distance,
        config::WEAK_HBOND.angle_degrees,
    ) {
        features |= FeatureBits::WEAK_HBOND;
    }

    features
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::contacts::find_contacts;
    use rspeggio_ccd::common::common_components;

    #[test]
    fn a_real_hydrophobic_contact_is_classified_between_two_ubq_side_chains() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");
        let components = common_components();
        let contacts = find_contacts(&pdb, config::CONTACT_TYPES_MAX_DIST);

        let found = contacts
            .iter()
            .any(|c| classify_features(c, &components).contains(FeatureBits::HYDROPHOBIC));
        assert!(
            found,
            "expected at least one hydrophobic contact in 1UBQ's core"
        );
    }

    #[test]
    fn a_real_carbonyl_contact_is_classified_somewhere_in_ubq() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");
        let components = common_components();
        let contacts = find_contacts(&pdb, config::CONTACT_TYPES_MAX_DIST);

        let found = contacts
            .iter()
            .any(|c| classify_features(c, &components).contains(FeatureBits::CARBONYL));
        assert!(
            found,
            "expected at least one backbone carbonyl-carbonyl contact in 1UBQ"
        );
    }

    #[test]
    fn a_contact_with_no_typed_atoms_has_no_features() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");
        let components: HashMap<String, CcdComponent> = HashMap::new();
        let contacts = find_contacts(&pdb, config::CONTACT_TYPES_MAX_DIST);
        assert!(
            !contacts.is_empty(),
            "sanity check: there are real contacts to classify"
        );

        for c in &contacts {
            assert!(classify_features(c, &components).is_empty());
        }
    }

    #[test]
    fn every_classified_contact_actually_respects_its_own_threshold() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");
        let components = common_components();
        let contacts = find_contacts(&pdb, config::CONTACT_TYPES_MAX_DIST);

        for c in &contacts {
            let features = classify_features(c, &components);
            if features.contains(FeatureBits::IONIC) {
                assert!(c.distance <= config::IONIC.distance);
            }
            if features.contains(FeatureBits::HYDROPHOBIC) {
                assert!(c.distance <= config::HYDROPHOBIC.distance);
            }
            if features.contains(FeatureBits::CARBONYL) {
                assert!(c.distance <= config::CARBONYL.distance);
            }
            if features.contains(FeatureBits::METAL_COMPLEX) {
                assert!(c.distance <= config::METAL.distance);
            }
            if features.contains(FeatureBits::HBOND) {
                assert!(c.distance <= config::HBOND.distance);
            }
            if features.contains(FeatureBits::WEAK_HBOND) {
                assert!(c.distance <= config::WEAK_HBOND.distance);
            }
            if features.contains(FeatureBits::AROMATIC) {
                assert!(c.distance <= config::AROMATIC.distance);
            }
        }
    }

    #[test]
    fn a_real_aromatic_contact_is_found_in_1ca2() {
        // 1UBQ's aromatic side chains never happen to pack within the 4.0A
        // threshold (confirmed directly: closest inter-residue
        // aromatic-aromatic atom pair is ~5.6A) -- 1CA2 is larger, with 238
        // aromatic ring atoms, and has real pairs as close as ~3.1A.
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1CA2.cif").expect("1CA2 should load");
        let components = common_components();
        let contacts = find_contacts(&pdb, config::CONTACT_TYPES_MAX_DIST);

        let found = contacts
            .iter()
            .any(|c| classify_features(c, &components).contains(FeatureBits::AROMATIC));
        assert!(found, "expected at least one real aromatic contact in 1CA2");
    }

    #[test]
    fn a_real_zinc_coordination_is_classified_as_metal_complex() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1CA2.cif").expect("1CA2 should load");
        let mut components = common_components();
        let zn =
            rspeggio_ccd::parser::load_ccd_component("../rspeggio-ccd/tests/fixtures/ZN_ideal.cif")
                .expect("ZN fixture should parse");
        components.insert("ZN".to_string(), zn);
        let contacts = find_contacts(&pdb, config::CONTACT_TYPES_MAX_DIST);

        let found = contacts
            .iter()
            .any(|c| classify_features(c, &components).contains(FeatureBits::METAL_COMPLEX));
        assert!(
            found,
            "expected at least one metal-complex contact around 1CA2's zinc"
        );
    }

    #[test]
    fn a_real_ionic_contact_is_found_somewhere_in_1ca2() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1CA2.cif").expect("1CA2 should load");
        let components = common_components();
        let contacts = find_contacts(&pdb, config::CONTACT_TYPES_MAX_DIST);

        let found = contacts
            .iter()
            .any(|c| classify_features(c, &components).contains(FeatureBits::IONIC));
        assert!(
            found,
            "expected at least one ionic (salt bridge) contact in 1CA2"
        );
    }

    #[test]
    fn a_real_hbond_is_found_via_placed_donor_hydrogen_geometry() {
        // BPTI (5PTI) is small, real, and has plenty of side-chain
        // donor/acceptor pairs (Ser/Thr/Asn/Gln/Lys/Arg/His side chains)
        // that don't need cross-residue lookup -- exactly the shapes
        // `donor_hydrogen_positions` covers.
        let (pdb, _errors) =
            pdbtbx::open("tests/fixtures/structures/5PTI.cif").expect("BPTI should load");
        let components = common_components();
        let contacts = find_contacts(&pdb, config::CONTACT_TYPES_MAX_DIST);

        let found = contacts
            .iter()
            .any(|c| classify_features(c, &components).contains(FeatureBits::HBOND));
        assert!(found, "expected at least one real hbond in BPTI");
    }

    #[test]
    fn a_real_weak_hbond_is_found_in_ubq() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");
        let components = common_components();
        let contacts = find_contacts(&pdb, config::CONTACT_TYPES_MAX_DIST);

        let found = contacts
            .iter()
            .any(|c| classify_features(c, &components).contains(FeatureBits::WEAK_HBOND));
        assert!(found, "expected at least one weak hbond (C-H...X) in 1UBQ");
    }
}
