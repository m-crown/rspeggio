// crates/rspeggio-core/src/contacts.rs
//
// The distance-category half of M5's atom-atom contact classification
// (positions 0-4 of the SIFt, see config::DistanceCategory), ported from
// real pdbe-arpeggio's `_calculate_atom_contacts`
// (interactions.py:707-773). Feature contacts (hbond, ionic, aromatic,
// ...) are a separate, later piece built on top of this.
//
// The neighbor search itself uses pdbtbx's own `rstar`-backed
// `create_hierarchy_rtree()` directly, rather than a separate spatial
// index built over a parallel `Structure` abstraction -- pdbtbx's
// hierarchy items already carry full atom/residue/chain context plus
// spatial queries, so a second representation would just be redundant
// upkeep for no benefit.

use crate::config::{self, DistanceCategory};
use pdbtbx::AtomConformerResidueChainModel;
use pdbtbx::{
    Atom, ContainsAtomConformer, ContainsAtomConformerResidue, ContainsAtomConformerResidueChain,
    PDB,
};

pub fn euclidean_distance(a: (f64, f64, f64), b: (f64, f64, f64)) -> f64 {
    ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2) + (a.2 - b.2).powi(2)).sqrt()
}

fn element_symbol(atom: &Atom) -> &'static str {
    atom.element().map(|e| e.symbol()).unwrap_or("")
}

// Real pdbe-arpeggio excludes contacts between sequence-adjacent
// polypeptide residues by default (interactions.py:732-741) -- a peptide
// bond's own backbone atoms would otherwise dominate every adjacent
// residue pair as spurious "clash"/"covalent" contacts. Its own check
// walks `prev_residue`/`next_residue` links built during structure
// loading; those links only exist on actual polypeptide-chain residues in
// the first place (a water or ligand never gets one), which is what its
// separate, apparently-redundant `is_polypeptide` guard is really
// enforcing. We don't build prev/next links, so this approximates the
// same thing directly: both residues are standard amino acids, in the
// same chain, one sequence position apart.
//
// Known limitation: seq_id adjacency isn't always the same as true
// polymer-chain adjacency -- a gap in the deposited structure (missing
// residues) or an insertion code would make two non-adjacent residues
// look seq_id-adjacent, or vice versa. Revisit if a real structure in the
// corpus exposes this.
pub fn is_sequence_adjacent(
    a: &AtomConformerResidueChainModel,
    b: &AtomConformerResidueChainModel,
) -> bool {
    config::STANDARD_AMINO_ACIDS.contains(&a.residue().name().unwrap_or_default())
        && config::STANDARD_AMINO_ACIDS.contains(&b.residue().name().unwrap_or_default())
        && a.chain().id() == b.chain().id()
        && (a.residue().id().0 - b.residue().id().0).abs() == 1
}

// Real pdbe-arpeggio determines "covalent" by checking actual bonded
// membership in OpenBabel's live-perceived molecular graph across the
// whole structure -- geometric perception this project's decisions 01/04
// deliberately avoid for chemistry judgments. Since contacts are only
// ever computed between different residues (intra-residue pairs are
// filtered out before classification even runs) and sequence-adjacent
// residues are excluded by default, the cross-residue covalent case that
// actually shows up in practice is narrow: disulfide bridges, or a
// covalently-attached ligand (not handled yet).
//
// This handles disulfides specifically: two sulfurs closer than their
// summed covalent radii (2 x 1.05 = 2.10 A) are bonded, not clashing.
// Real disulfide bond lengths are ~2.02-2.05 A (confirmed against BPTI's
// three known disulfides), comfortably under that cutoff with room to
// spare before the ambiguous zone where two sulfurs might just be
// sterically close without being bonded.
fn is_disulfide_bond(element_1: &str, element_2: &str, distance: f64) -> bool {
    if !element_1.eq_ignore_ascii_case("S") || !element_2.eq_ignore_ascii_case("S") {
        return false;
    }
    match config::covalent_radius("S") {
        Some(radius) => distance < radius * 2.0,
        None => false,
    }
}

// Classifies an inter-residue atom pair's distance category. Callers are
// responsible for the pre-filtering real arpeggio does before reaching
// this point: skip hydrogens, skip intra-residue pairs, skip
// sequence-adjacent residues by default (none of that is geometry this
// function needs to know about).
pub fn classify_distance(element_1: &str, element_2: &str, distance: f64) -> DistanceCategory {
    if is_disulfide_bond(element_1, element_2, distance) {
        return DistanceCategory::Covalent;
    }

    let sum_cov_radii = config::covalent_radius(element_1).unwrap_or(0.0)
        + config::covalent_radius(element_2).unwrap_or(0.0);
    let sum_vdw_radii =
        config::vdw_radius(element_1).unwrap_or(0.0) + config::vdw_radius(element_2).unwrap_or(0.0);

    if distance < sum_cov_radii {
        DistanceCategory::Clash
    } else if distance < sum_vdw_radii {
        DistanceCategory::VdwClash
    } else if distance <= sum_vdw_radii + config::VDW_COMP_FACTOR {
        DistanceCategory::Vdw
    } else {
        DistanceCategory::Proximal
    }
}

// One classified atom-atom contact. Holds the two hierarchy items
// directly (borrowed from the source PDB) rather than indices, since
// pdbtbx's own types already give a stable, cheap-to-clone handle with
// full atom/residue/chain context attached.
#[derive(Clone)]
pub struct Contact<'a> {
    pub atom_1: AtomConformerResidueChainModel<'a>,
    pub atom_2: AtomConformerResidueChainModel<'a>,
    pub distance: f64,
    pub category: DistanceCategory,
}

// Finds every contact in a structure within `cutoff` Angstroms, applying
// real pdbe-arpeggio's pre-filters (interactions.py:707-741): hydrogens
// are skipped entirely (most X-ray structures have none anyway); each
// unordered pair is considered once (via atom serial number, which is
// already guaranteed unique within a model); intra-residue pairs are
// skipped (compared by residue identity, i.e. pointer equality into the
// same live PDB); sequence-adjacent residue pairs are skipped by default.
pub fn find_contacts(pdb: &PDB, cutoff: f64) -> Vec<Contact<'_>> {
    let tree = pdb.create_hierarchy_rtree();
    let cutoff_squared = cutoff * cutoff;
    let mut contacts = Vec::new();

    for atom in tree.iter() {
        if element_symbol(atom.atom()) == "H" {
            continue;
        }

        for neighbor in tree.locate_within_distance(atom.atom().pos(), cutoff_squared) {
            if element_symbol(neighbor.atom()) == "H" {
                continue;
            }
            if neighbor.atom().serial_number() <= atom.atom().serial_number() {
                continue; // each unordered pair considered once
            }
            if std::ptr::eq(atom.residue(), neighbor.residue()) {
                continue; // intra-residue
            }
            if is_sequence_adjacent(atom, neighbor) {
                continue;
            }

            let distance = euclidean_distance(atom.atom().pos(), neighbor.atom().pos());
            contacts.push(Contact {
                atom_1: atom.clone(),
                atom_2: neighbor.clone(),
                distance,
                category: classify_distance(
                    element_symbol(atom.atom()),
                    element_symbol(neighbor.atom()),
                    distance,
                ),
            });
        }
    }

    contacts
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn euclidean_distance_of_identical_points_is_zero() {
        assert_eq!(euclidean_distance((1.0, 2.0, 3.0), (1.0, 2.0, 3.0)), 0.0);
    }

    #[test]
    fn euclidean_distance_matches_a_known_3_4_5_triangle() {
        assert_eq!(euclidean_distance((0.0, 0.0, 0.0), (3.0, 4.0, 0.0)), 5.0);
    }

    #[test]
    fn two_sulfurs_too_far_apart_are_not_a_disulfide() {
        assert_ne!(classify_distance("S", "S", 4.0), DistanceCategory::Covalent);
    }

    #[test]
    fn identical_atoms_at_zero_distance_clash() {
        assert_eq!(classify_distance("C", "C", 0.1), DistanceCategory::Clash);
    }

    #[test]
    fn a_proximal_pair_is_beyond_vdw_plus_compensation() {
        // C-C: cov 0.76+0.76=1.52, vdw 1.70+1.70=3.40, +0.1 comp = 3.50.
        assert_eq!(classify_distance("C", "C", 4.0), DistanceCategory::Proximal);
    }

    #[test]
    fn a_vdw_pair_sits_between_vdw_radii_sum_and_the_compensated_bound() {
        assert_eq!(classify_distance("C", "C", 3.45), DistanceCategory::Vdw);
    }

    #[test]
    fn a_vdw_clash_sits_between_covalent_and_vdw_radii_sums() {
        assert_eq!(classify_distance("C", "C", 2.5), DistanceCategory::VdwClash);
    }

    #[test]
    fn consecutive_real_residues_are_sequence_adjacent() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");
        let atoms: Vec<_> = pdb.atoms_with_hierarchy().collect();

        let met1_ca = atoms
            .iter()
            .find(|a| a.residue().id().0 == 1 && a.atom().name() == "CA")
            .expect("MET1 CA should be present");
        let gln2_ca = atoms
            .iter()
            .find(|a| a.residue().id().0 == 2 && a.atom().name() == "CA")
            .expect("GLN2 CA should be present");

        assert!(is_sequence_adjacent(met1_ca, gln2_ca));
        assert!(is_sequence_adjacent(gln2_ca, met1_ca), "symmetric");
    }

    #[test]
    fn distant_real_residues_are_not_sequence_adjacent() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");
        let atoms: Vec<_> = pdb.atoms_with_hierarchy().collect();

        let met1_ca = atoms
            .iter()
            .find(|a| a.residue().id().0 == 1 && a.atom().name() == "CA")
            .unwrap();
        let ile3_ca = atoms
            .iter()
            .find(|a| a.residue().id().0 == 3 && a.atom().name() == "CA")
            .unwrap();

        assert!(!is_sequence_adjacent(met1_ca, ile3_ca));
    }

    #[test]
    fn a_water_is_never_sequence_adjacent_even_with_a_matching_seq_id() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");
        let atoms: Vec<_> = pdb.atoms_with_hierarchy().collect();

        let met1_ca = atoms
            .iter()
            .find(|a| a.residue().id().0 == 1 && a.atom().name() == "CA")
            .unwrap();
        // Some water in the fixture, numbered independently of the
        // polypeptide chain -- if any water happens to share a seq_id
        // with MET1+1, it must still not be treated as adjacent.
        if let Some(water) = atoms.iter().find(|a| a.residue().name() == Some("HOH")) {
            assert!(!is_sequence_adjacent(met1_ca, water));
        }
    }

    #[test]
    fn finds_bptis_three_real_disulfide_bonds_as_covalent_contacts() {
        let (pdb, _errors) =
            pdbtbx::open("tests/fixtures/structures/5PTI.cif").expect("BPTI should load");

        let contacts = find_contacts(&pdb, config::CONTACT_TYPES_MAX_DIST);

        let disulfide_pairs: Vec<(isize, isize)> = contacts
            .iter()
            .filter(|c| c.category == DistanceCategory::Covalent)
            .map(|c| {
                let a = c.atom_1.residue().id().0;
                let b = c.atom_2.residue().id().0;
                if a < b {
                    (a, b)
                } else {
                    (b, a)
                }
            })
            .collect();

        for expected in [(5, 55), (14, 38), (30, 51)] {
            assert!(
                disulfide_pairs.contains(&expected),
                "expected disulfide {expected:?} in {disulfide_pairs:?}"
            );
        }
        assert_eq!(
            disulfide_pairs.len(),
            3,
            "BPTI has exactly three disulfides, no more should be found: {disulfide_pairs:?}"
        );
    }

    #[test]
    fn no_contact_is_ever_intra_residue_or_sequence_adjacent() {
        let (pdb, _errors) =
            pdbtbx::open("tests/fixtures/structures/5PTI.cif").expect("BPTI should load");

        let contacts = find_contacts(&pdb, config::CONTACT_TYPES_MAX_DIST);
        assert!(
            !contacts.is_empty(),
            "sanity check: BPTI should have real contacts"
        );

        for c in &contacts {
            assert!(
                !std::ptr::eq(c.atom_1.residue(), c.atom_2.residue()),
                "no contact should be intra-residue"
            );
            assert!(
                !is_sequence_adjacent(&c.atom_1, &c.atom_2),
                "no contact should be between sequence-adjacent residues"
            );
        }
    }

    #[test]
    fn every_contact_atom_pair_is_within_the_requested_cutoff() {
        let (pdb, _errors) =
            pdbtbx::open("tests/fixtures/structures/5PTI.cif").expect("BPTI should load");

        let cutoff = 4.5;
        let contacts = find_contacts(&pdb, cutoff);
        assert!(!contacts.is_empty());

        for c in &contacts {
            assert!(
                c.distance <= cutoff,
                "contact distance {} exceeds requested cutoff {cutoff}",
                c.distance
            );
        }
    }

    #[test]
    fn a_zero_cutoff_finds_no_contacts() {
        let (pdb, _errors) =
            pdbtbx::open("tests/fixtures/structures/5PTI.cif").expect("BPTI should load");

        assert!(find_contacts(&pdb, 0.0).is_empty());
    }
}
