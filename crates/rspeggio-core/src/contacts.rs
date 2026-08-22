// crates/rspeggio-core/src/contacts.rs
//
// The distance-category half of M5's atom-atom contact classification
// (positions 0-4 of the SIFt, see config::DistanceCategory), ported from
// real pdbe-arpeggio's `_calculate_atom_contacts`
// (interactions.py:707-773). Feature contacts (hbond, ionic, aromatic,
// ...) are a separate, later piece built on top of this.

use crate::config::{self, DistanceCategory};
use crate::structure::StructureResidue;

pub fn euclidean_distance(a: (f64, f64, f64), b: (f64, f64, f64)) -> f64 {
    ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2) + (a.2 - b.2).powi(2)).sqrt()
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
pub fn is_sequence_adjacent(a: &StructureResidue, b: &StructureResidue) -> bool {
    config::STANDARD_AMINO_ACIDS.contains(&a.comp_id.as_str())
        && config::STANDARD_AMINO_ACIDS.contains(&b.comp_id.as_str())
        && a.chain_id == b.chain_id
        && (a.seq_id - b.seq_id).abs() == 1
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
// three known disulfides, see typing tests), comfortably under that
// cutoff with room to spare before the ambiguous zone where two sulfurs
// might just be sterically close without being bonded.
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::load_structure;

    #[test]
    fn euclidean_distance_of_identical_points_is_zero() {
        assert_eq!(euclidean_distance((1.0, 2.0, 3.0), (1.0, 2.0, 3.0)), 0.0);
    }

    fn ubq_residue(structure: &crate::structure::Structure, seq_id: isize) -> StructureResidue {
        // Cloning by hand since StructureResidue doesn't derive Clone --
        // fine for tests, avoids adding a Clone impl the real pipeline
        // never needs.
        let r = structure
            .residues
            .iter()
            .find(|r| r.seq_id == seq_id)
            .unwrap_or_else(|| panic!("residue {seq_id} should be present"));
        StructureResidue {
            comp_id: r.comp_id.clone(),
            chain_id: r.chain_id.clone(),
            seq_id: r.seq_id,
            atoms: Vec::new(),
        }
    }

    #[test]
    fn consecutive_real_residues_are_sequence_adjacent() {
        let structure =
            load_structure("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");

        // MET1, GLN2 -- confirmed real consecutive residues in the fixture.
        let met1 = ubq_residue(&structure, 1);
        let gln2 = ubq_residue(&structure, 2);

        assert!(is_sequence_adjacent(&met1, &gln2));
        assert!(is_sequence_adjacent(&gln2, &met1), "symmetric");
    }

    #[test]
    fn distant_real_residues_are_not_sequence_adjacent() {
        let structure =
            load_structure("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");

        let met1 = ubq_residue(&structure, 1);
        let ile3 = ubq_residue(&structure, 3);

        assert!(!is_sequence_adjacent(&met1, &ile3));
    }

    #[test]
    fn a_water_is_never_sequence_adjacent_even_with_a_matching_seq_id() {
        let structure =
            load_structure("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");

        let met1 = ubq_residue(&structure, 1);
        let water = StructureResidue {
            comp_id: "HOH".to_string(),
            chain_id: met1.chain_id.clone(),
            seq_id: 2, // numerically adjacent to MET1, but not a polypeptide residue
            atoms: Vec::new(),
        };

        assert!(!is_sequence_adjacent(&met1, &water));
    }

    #[test]
    fn residues_in_different_chains_are_never_sequence_adjacent() {
        let a = StructureResidue {
            comp_id: "ALA".to_string(),
            chain_id: "A".to_string(),
            seq_id: 1,
            atoms: Vec::new(),
        };
        let b = StructureResidue {
            comp_id: "GLY".to_string(),
            chain_id: "B".to_string(),
            seq_id: 2,
            atoms: Vec::new(),
        };

        assert!(!is_sequence_adjacent(&a, &b));
    }

    #[test]
    fn euclidean_distance_matches_a_known_3_4_5_triangle() {
        assert_eq!(euclidean_distance((0.0, 0.0, 0.0), (3.0, 4.0, 0.0)), 5.0);
    }

    #[test]
    fn bptis_real_disulfide_bonds_classify_as_covalent() {
        let structure =
            load_structure("tests/fixtures/structures/5PTI.cif").expect("BPTI fixture should load");

        let sg = |seq_id: isize| -> (f64, f64, f64) {
            structure
                .residues
                .iter()
                .find(|r| r.comp_id == "CYS" && r.seq_id == seq_id)
                .and_then(|r| r.atoms.iter().find(|a| a.name == "SG"))
                .unwrap_or_else(|| panic!("CYS{seq_id} SG should be present"))
                .pos
        };

        // BPTI's three known disulfides: 5-55, 14-38, 30-51.
        for (a, b) in [(5, 55), (14, 38), (30, 51)] {
            let distance = euclidean_distance(sg(a), sg(b));
            assert_eq!(
                classify_distance("S", "S", distance),
                DistanceCategory::Covalent,
                "CYS{a}-CYS{b} SG-SG at {distance:.3} A should be a disulfide bond"
            );
        }
    }

    #[test]
    fn two_sulfurs_too_far_apart_are_not_a_disulfide() {
        // Same element pair as a real disulfide, but a distance
        // (arbitrarily) beyond bonding range -- should fall through to
        // ordinary vdw-radii-based classification instead.
        assert_ne!(classify_distance("S", "S", 4.0), DistanceCategory::Covalent);
    }

    #[test]
    fn identical_atoms_at_zero_distance_clash() {
        // Not a realistic contact (would be filtered as intra-residue or
        // even the same atom upstream), but confirms the boundary logic:
        // well under any covalent-radii sum falls to Clash for a
        // non-disulfide pair.
        assert_eq!(classify_distance("C", "C", 0.1), DistanceCategory::Clash);
    }

    #[test]
    fn a_proximal_pair_is_beyond_vdw_plus_compensation() {
        // C-C: cov 0.76+0.76=1.52, vdw 1.70+1.70=3.40, +0.1 comp = 3.50.
        assert_eq!(classify_distance("C", "C", 4.0), DistanceCategory::Proximal);
    }

    #[test]
    fn a_vdw_pair_sits_between_vdw_radii_sum_and_the_compensated_bound() {
        // Between 3.40 (vdw sum) and 3.50 (vdw sum + compensation).
        assert_eq!(classify_distance("C", "C", 3.45), DistanceCategory::Vdw);
    }

    #[test]
    fn a_vdw_clash_sits_between_covalent_and_vdw_radii_sums() {
        // Between 1.52 (cov sum) and 3.40 (vdw sum).
        assert_eq!(classify_distance("C", "C", 2.5), DistanceCategory::VdwClash);
    }
}
