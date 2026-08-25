// crates/rspeggio-core/src/unsatisfied.rs
//
// M7b: unsatisfied-contact detection, real pdbe-arpeggio's per-atom
// "potential vs actual" SIFt matching (`_initialize_atom_sift`'s
// `potential_fsift`, `update_atom_fsift`, and `utils.sift_match_base3`/
// `human_sift_match`). Given the atom-atom feature contacts M5 already
// computes, this answers a different question: for one real atom, out of
// every feature its own chemistry makes it *capable* of (an acceptor with
// no partner, a hydrophobe buried with nothing nearby to pack against,
// ...), which ones did a real contact actually satisfy, and which never
// found a partner at all? "Unsatisfied" contacts (an atom capable of
// something but never realizing it) are exactly the ones with no
// partner -- additive to, not replacing, the exact-parity JSON export
// (matches the original project plan's own framing of this milestone).
//
// Scope decision, following this project's established pattern of tight
// initial slices: whole-structure only, no selection scoping (unlike
// `export.rs`, which real arpeggio -- and this project -- restricts to
// `selection_plus`). Real arpeggio's own per-atom SIFt matching iterates
// `self.selection_plus` too, so this is a real, deliberate simplification,
// not a port of the oracle's exact scope -- worth adding selection
// awareness later if it's needed, not attempted here. Also not ported:
// real arpeggio's numeric `potential_hbonds`/`potential_polars` lone-pair
// counts (`config.VALENCE` arithmetic) -- this only reproduces the binary
// per-feature capability check `sift_match_base3` actually needs.

use crate::config::{self, FeatureBits};
use crate::contacts::find_contacts;
use crate::export::FEATURE_LABELS;
use crate::features::classify_features;
use crate::join::join_atom;
use crate::typing::{type_atom, AtomTypeBits};
use pdbtbx::{
    Atom, AtomConformerResidueChainModel, ContainsAtomConformer, ContainsAtomConformerResidue, PDB,
};
use rspeggio_ccd::component::CcdComponent;
use std::collections::HashMap;

// Real pdbe-arpeggio's per-atom "potential" feature capacity
// (`_initialize_atom_sift`'s `potential_fsift`, `interactions.py:1791-1852`):
// whether this atom's own typing makes it *capable* of a given feature at
// all, independent of whether any real contact actually realizes it. Uses
// the same bit positions/order as `FeatureBits` (this project's SIFt
// feature order already matches the oracle's).
pub fn potential_features(bits: AtomTypeBits, element: &str) -> FeatureBits {
    let mut potential = FeatureBits::empty();

    let hbond_capable =
        bits.contains(AtomTypeBits::HBOND_ACCEPTOR) || bits.contains(AtomTypeBits::HBOND_DONOR);
    if hbond_capable {
        potential |= FeatureBits::HBOND | FeatureBits::POLAR;
    }

    // Real arpeggio's own `is_halogen` check (`config.HALOGENS`) is
    // broader than this project's `WEAK_HBOND_ACCEPTOR` bit alone (which
    // deliberately excludes a bare halide ion -- see `typing.rs`'s own
    // note on why): a halide ion still has real weak-hbond *potential* by
    // this definition, even though it isn't typed as an actual acceptor.
    let is_halogen = matches!(
        element.to_uppercase().as_str(),
        "F" | "CL" | "BR" | "I" | "AT"
    );
    let weak_hbond_capable = bits.contains(AtomTypeBits::WEAK_HBOND_ACCEPTOR)
        || bits.contains(AtomTypeBits::WEAK_HBOND_DONOR)
        || hbond_capable
        || is_halogen;
    if weak_hbond_capable {
        potential |= FeatureBits::WEAK_HBOND | FeatureBits::WEAK_POLAR;
    }

    // Real "xbond acceptor" is the same atom set as "hbond acceptor" (see
    // `typing.rs`'s own note on this).
    if bits.contains(AtomTypeBits::HBOND_ACCEPTOR) || bits.contains(AtomTypeBits::XBOND_DONOR) {
        potential |= FeatureBits::XBOND;
    }

    if bits.contains(AtomTypeBits::POS_IONISABLE) || bits.contains(AtomTypeBits::NEG_IONISABLE) {
        potential |= FeatureBits::IONIC;
    }

    if bits.contains(AtomTypeBits::HBOND_ACCEPTOR) || bits.contains(AtomTypeBits::METAL) {
        potential |= FeatureBits::METAL_COMPLEX;
    }

    if bits.contains(AtomTypeBits::AROMATIC) {
        potential |= FeatureBits::AROMATIC;
    }

    if bits.contains(AtomTypeBits::HYDROPHOBE) {
        potential |= FeatureBits::HYDROPHOBIC;
    }

    if bits.contains(AtomTypeBits::CARBONYL_OXYGEN) || bits.contains(AtomTypeBits::CARBONYL_CARBON)
    {
        potential |= FeatureBits::CARBONYL;
    }

    potential
}

// Real `utils.sift_match_base3`'s per-position outcome. Real arpeggio
// treats "not potential but actually observed" as an invariant violation
// (raises `SiftMatchError`) -- this instead just reports `Matched`
// (`actual` is the more authoritative signal if the two ever disagree,
// and a hard error here would make one inconsistent atom fail the whole
// computation rather than just this one entry).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiftMatch {
    Unmatched,   // capable of this feature, but no real contact ever satisfied it
    Matched,     // capable, and at least one real contact did satisfy it
    NotPossible, // this atom's own typing has no capability for this feature at all
}

fn match_feature(potential: FeatureBits, actual: FeatureBits, bit: FeatureBits) -> SiftMatch {
    if actual.contains(bit) {
        SiftMatch::Matched
    } else if potential.contains(bit) {
        SiftMatch::Unmatched
    } else {
        SiftMatch::NotPossible
    }
}

// One real atom's potential-vs-actual feature capability, across every
// real contact it participates in anywhere in the structure.
pub struct AtomSiftMatch<'a> {
    pub hierarchy: AtomConformerResidueChainModel<'a>,
    pub potential: FeatureBits,
    pub actual: FeatureBits,
}

impl AtomSiftMatch<'_> {
    pub fn match_code(&self, bit: FeatureBits) -> SiftMatch {
        match_feature(self.potential, self.actual, bit)
    }

    // Real arpeggio's `utils.human_sift_match`: e.g.
    // `"Matched hbond:Unmatched ionic"` -- skips any feature this atom has
    // no potential for at all, alphabetically sorted (matching the
    // oracle's own `terms.sort()`).
    pub fn human_readable(&self) -> String {
        let mut terms: Vec<String> = FEATURE_LABELS
            .iter()
            .filter_map(
                |(bit, label)| match match_feature(self.potential, self.actual, *bit) {
                    SiftMatch::NotPossible => None,
                    SiftMatch::Matched => Some(format!("Matched {label}")),
                    SiftMatch::Unmatched => Some(format!("Unmatched {label}")),
                },
            )
            .collect();
        terms.sort();
        terms.join(":")
    }
}

// Every real atom's potential-vs-actual feature match, across the whole
// structure. `Err` as soon as any residue's comp_id has no matching
// `CcdComponent` -- consistent with decision 03 (unknown components fail
// loudly): a real atom's `potential` capability genuinely can't be
// computed without its CCD typing.
pub fn compute_atom_sift_matches<'a>(
    pdb: &'a PDB,
    components: &'a HashMap<String, CcdComponent>,
) -> Result<Vec<AtomSiftMatch<'a>>, String> {
    let mut matches: HashMap<*const Atom, AtomSiftMatch<'a>> = HashMap::new();

    for hierarchy in pdb.atoms_with_hierarchy() {
        let joined = join_atom(hierarchy.clone(), components);
        let (Some(ccd_atom), Some(component)) = (joined.ccd_atom, joined.component) else {
            let comp_id = hierarchy.residue().name().unwrap_or_default();
            return Err(format!("unknown CCD component: {comp_id}"));
        };
        let bits = type_atom(ccd_atom, component);
        let potential = potential_features(bits, ccd_atom.element());
        matches.insert(
            hierarchy.atom() as *const Atom,
            AtomSiftMatch {
                hierarchy,
                potential,
                actual: FeatureBits::empty(),
            },
        );
    }

    for contact in find_contacts(pdb, config::CONTACT_TYPES_MAX_DIST) {
        let features = classify_features(&contact, components);
        if features.is_empty() {
            continue;
        }
        if let Some(entry) = matches.get_mut(&(contact.atom_1.atom() as *const Atom)) {
            entry.actual |= features;
        }
        if let Some(entry) = matches.get_mut(&(contact.atom_2.atom() as *const Atom)) {
            entry.actual |= features;
        }
    }

    let mut result: Vec<AtomSiftMatch<'a>> = matches.into_values().collect();
    result.sort_by_key(|m| m.hierarchy.atom().serial_number());
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use pdbtbx::ContainsAtomConformerResidue;
    use rspeggio_ccd::common::common_components;

    fn find_match<'a>(
        matches: &'a [AtomSiftMatch<'a>],
        comp_id: &str,
        seq_id: isize,
        atom_name: &str,
    ) -> &'a AtomSiftMatch<'a> {
        matches
            .iter()
            .find(|m| {
                m.hierarchy.residue().name() == Some(comp_id)
                    && m.hierarchy.residue().id().0 == seq_id
                    && m.hierarchy.atom().name() == atom_name
            })
            .unwrap_or_else(|| panic!("{comp_id}{seq_id} {atom_name} should be present"))
    }

    #[test]
    fn a_real_zinc_coordinating_histidine_nitrogen_has_a_matched_metal_complex() {
        // HIS94's real, specific zinc-coordinating atom is NE2 (confirmed
        // directly via an ad-hoc sweep of `classify_features`'
        // METAL_COMPLEX contacts in 1CA2 before writing this test, at
        // ~1.99A -- not ND1, even though ND1 is also within the broader
        // 4.5A atom-atom cutoff for *some* classification). Its own
        // METAL_COMPLEX potential should be genuinely satisfied here.
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1CA2.cif").expect("1CA2 should load");
        let mut components = common_components();
        let zn =
            rspeggio_ccd::parser::load_ccd_component("../rspeggio-ccd/tests/fixtures/ZN_ideal.cif")
                .expect("ZN fixture should parse");
        components.insert("ZN".to_string(), zn);

        let matches = compute_atom_sift_matches(&pdb, &components)
            .expect("1CA2 is amino acids + water + zinc, all now known");

        let his_ne2 = find_match(&matches, "HIS", 94, "NE2");
        assert_eq!(
            his_ne2.match_code(FeatureBits::METAL_COMPLEX),
            SiftMatch::Matched,
            "HIS94 NE2 really coordinates 1CA2's zinc"
        );
    }

    #[test]
    fn a_plain_hydrophobic_carbon_has_no_ionic_potential_at_all() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");
        let components = common_components();

        let matches =
            compute_atom_sift_matches(&pdb, &components).expect("1UBQ should all be known");

        // ALA's CB is a plain methyl carbon -- real chemistry, no ionisable
        // capacity whatsoever (confirmed already by `typing.rs`'s own
        // tests for the same atom).
        let ala1_cb = matches
            .iter()
            .find(|m| {
                m.hierarchy.residue().name() == Some("ALA") && m.hierarchy.atom().name() == "CB"
            })
            .expect("1UBQ should have at least one ALA CB");
        assert_eq!(
            ala1_cb.match_code(FeatureBits::IONIC),
            SiftMatch::NotPossible
        );
    }

    #[test]
    fn a_real_unsatisfied_hbond_acceptor_exists_somewhere_in_1ubq() {
        // Confirmed via an ad-hoc sweep before writing this test: not
        // every real hbond-acceptor-capable atom in a real protein finds
        // a real partner within cutoff -- surface-exposed carbonyls
        // pointing into open solvent are the classic real case.
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");
        let components = common_components();

        let matches =
            compute_atom_sift_matches(&pdb, &components).expect("1UBQ should all be known");

        let found = matches
            .iter()
            .any(|m| m.match_code(FeatureBits::HBOND) == SiftMatch::Unmatched);
        assert!(
            found,
            "expected at least one real hbond-capable atom with no satisfied hbond in 1UBQ"
        );
    }

    #[test]
    fn human_readable_output_matches_real_arpeggios_format() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1CA2.cif").expect("1CA2 should load");
        let mut components = common_components();
        let zn =
            rspeggio_ccd::parser::load_ccd_component("../rspeggio-ccd/tests/fixtures/ZN_ideal.cif")
                .expect("ZN fixture should parse");
        components.insert("ZN".to_string(), zn);

        let matches = compute_atom_sift_matches(&pdb, &components).expect("1CA2 should export");
        let his_ne2 = find_match(&matches, "HIS", 94, "NE2");

        let human = his_ne2.human_readable();
        assert!(
            human.contains("Matched metal_complex"),
            "expected 'Matched metal_complex' in {human:?}"
        );
        // No entry for a feature this atom has zero potential for
        // (e.g. carbonyl -- NE2 is a ring nitrogen, not a carbonyl atom).
        assert!(!human.contains("carbonyl"));
    }

    #[test]
    fn an_unknown_component_fails_the_whole_computation_loudly() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");
        let components: HashMap<String, CcdComponent> = HashMap::new();

        assert!(compute_atom_sift_matches(&pdb, &components).is_err());
    }
}
