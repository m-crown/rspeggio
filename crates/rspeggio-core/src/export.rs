// crates/rspeggio-core/src/export.rs
//
// M7: JSON export matching real pdbe-arpeggio's own schema
// (`Arpeggio.get_contacts` in `interactions.py`), for the atom-atom
// contact type only -- the plane-plane/atom-plane/group-plane types
// `rings.rs` already computes the geometry for aren't wired into export
// yet, a follow-up slice, not a silent gap.
//
// Scope decision (asked, not guessed): this only implements
// "whole-structure" mode -- every residue treated as if it were real
// pdbe-arpeggio's user-supplied `-s` selection. Real arpeggio's
// `interacting_entities` field (INTER/INTRA_SELECTION/INTRA_NON_SELECTION/
// SELECTION_WATER/NON_SELECTION_WATER/WATER_WATER) depends entirely on
// that selection, which this project has no concept of at all yet.
// Working out what each `tests/fixtures/golden/*.json` file's real
// selection actually was (to reproduce it exactly) is real, separate scope
// -- deferred, not attempted here. With everything selected, real
// arpeggio's own `__get_contact_type` logic collapses to exactly three
// outcomes (verified by hand-tracing its real if-statement order, not
// guessed): both real atoms water -> WATER_WATER; exactly one water ->
// SELECTION_WATER; neither -> INTRA_SELECTION. INTER/INTRA_NON_SELECTION/
// NON_SELECTION_WATER never occur in this mode (they all require an atom
// *outside* the selection, which doesn't exist here).
//
// This does NOT byte-match the golden fixtures' own `interacting_entities`
// values (those used a real, different selection) -- but `bgn`/`end`/
// `distance`/`contact` for a given real atom pair are independent of
// selection, so those fields *are* checkable against golden data directly,
// and this module's tests do exactly that.

use crate::config::{self, DistanceCategory, FeatureBits};
use crate::contacts::{find_contacts, Contact};
use crate::features::classify_features;
use pdbtbx::{
    AtomConformerResidueChainModel, ContainsAtomConformer, ContainsAtomConformerResidue,
    ContainsAtomConformerResidueChain, PDB,
};
use rspeggio_ccd::component::{CcdComponent, ComponentType};
use serde::Serialize;
use std::collections::HashMap;

#[derive(Debug, Clone, Serialize)]
pub struct AtomIdentity {
    pub auth_asym_id: String,
    pub auth_atom_id: String,
    pub auth_seq_id: isize,
    pub label_comp_id: String,
    pub label_comp_type: &'static str,
    #[serde(rename = "pdbx_PDB_ins_code")]
    pub pdbx_pdb_ins_code: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct AtomAtomContactJson {
    pub bgn: AtomIdentity,
    pub end: AtomIdentity,
    #[serde(rename = "type")]
    pub entry_type: &'static str,
    pub distance: f64,
    pub contact: Vec<&'static str>,
    pub interacting_entities: &'static str,
}

// Real pdbe-arpeggio's own 15-position label list
// (`Arpeggio.get_contacts`, `interactions.py:178`) -- the first 5 are
// `DistanceCategory`, mutually exclusive (exactly one always present); the
// remaining 10 are `FeatureBits`, zero or more present. The feature order
// here matches `FeatureBits`' own declared bit order exactly (both are the
// same real SIFt position order), so this list doubles as a decode table
// for the bitflags themselves.
const FEATURE_LABELS: &[(FeatureBits, &str)] = &[
    (FeatureBits::HBOND, "hbond"),
    (FeatureBits::WEAK_HBOND, "weak_hbond"),
    (FeatureBits::XBOND, "xbond"),
    (FeatureBits::IONIC, "ionic"),
    (FeatureBits::METAL_COMPLEX, "metal_complex"),
    (FeatureBits::AROMATIC, "aromatic"),
    (FeatureBits::HYDROPHOBIC, "hydrophobic"),
    (FeatureBits::CARBONYL, "carbonyl"),
    (FeatureBits::POLAR, "polar"),
    (FeatureBits::WEAK_POLAR, "weak_polar"),
];

fn distance_category_label(category: DistanceCategory) -> &'static str {
    match category {
        DistanceCategory::Clash => "clash",
        DistanceCategory::Covalent => "covalent",
        DistanceCategory::VdwClash => "vdw_clash",
        DistanceCategory::Vdw => "vdw",
        DistanceCategory::Proximal => "proximal",
    }
}

fn contact_labels(category: DistanceCategory, features: FeatureBits) -> Vec<&'static str> {
    let mut labels = vec![distance_category_label(category)];
    for (bit, label) in FEATURE_LABELS {
        if features.contains(*bit) {
            labels.push(label);
        }
    }
    labels
}

// Real pdbe-arpeggio rounds with `round(np.float64(distance), 2)`
// (`interactions.py:190`), Python/numpy's round-half-to-even. This uses
// ordinary round-half-away-from-zero instead (`f64::round`) -- the two
// only disagree on an exact tie at the second decimal place, vanishingly
// unlikely for a real, continuously-varying Euclidean distance, but a real
// (if practically unobservable) divergence worth naming rather than
// silently assuming away.
fn round_2(x: f64) -> f64 {
    (x * 100.0).round() / 100.0
}

fn insertion_code(residue: &pdbtbx::Residue) -> String {
    match residue.id().1 {
        Some(code) => code.to_string(),
        None => " ".to_string(),
    }
}

fn atom_identity(atom: &AtomConformerResidueChainModel, component: &CcdComponent) -> AtomIdentity {
    let residue = atom.residue();
    AtomIdentity {
        auth_asym_id: atom.chain().id().to_string(),
        auth_atom_id: atom.atom().name().to_string(),
        auth_seq_id: residue.id().0,
        label_comp_id: residue.name().unwrap_or_default().to_string(),
        label_comp_type: component.component_type().code(),
        pdbx_pdb_ins_code: insertion_code(residue),
    }
}

// Real `__get_contact_type`'s logic (`interactions.py:643`), specialized
// to "every residue is in the selection" -- see module doc for the
// hand-traced derivation of why only these 3 outcomes remain.
fn whole_structure_interacting_entities(a: ComponentType, b: ComponentType) -> &'static str {
    match (a == ComponentType::Water, b == ComponentType::Water) {
        (true, true) => "WATER_WATER",
        (true, false) | (false, true) => "SELECTION_WATER",
        (false, false) => "INTRA_SELECTION",
    }
}

fn component_for<'a>(
    atom: &AtomConformerResidueChainModel,
    components: &'a HashMap<String, CcdComponent>,
) -> Result<&'a CcdComponent, String> {
    let comp_id = atom.residue().name().unwrap_or_default();
    components
        .get(comp_id)
        .ok_or_else(|| format!("unknown CCD component: {comp_id}"))
}

fn build_atom_atom_entry(
    contact: &Contact,
    components: &HashMap<String, CcdComponent>,
) -> Result<AtomAtomContactJson, String> {
    let component_1 = component_for(&contact.atom_1, components)?;
    let component_2 = component_for(&contact.atom_2, components)?;

    let features = classify_features(contact, components);
    let interacting_entities = whole_structure_interacting_entities(
        component_1.component_type(),
        component_2.component_type(),
    );

    Ok(AtomAtomContactJson {
        bgn: atom_identity(&contact.atom_1, component_1),
        end: atom_identity(&contact.atom_2, component_2),
        entry_type: "atom-atom",
        distance: round_2(contact.distance),
        contact: contact_labels(contact.category, features),
        interacting_entities,
    })
}

// Every atom-atom contact in `pdb`, in real pdbe-arpeggio's JSON shape,
// whole-structure mode (see module doc). `Err` as soon as any residue's
// comp_id has no matching `CcdComponent` -- consistent with decision 03
// (unknown components fail loudly), since `label_comp_type` genuinely
// can't be produced without one.
pub fn export_atom_atom_contacts(
    pdb: &PDB,
    components: &HashMap<String, CcdComponent>,
) -> Result<Vec<AtomAtomContactJson>, String> {
    find_contacts(pdb, config::CONTACT_TYPES_MAX_DIST)
        .iter()
        .map(|c| build_atom_atom_entry(c, components))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rspeggio_ccd::common::common_components;

    #[test]
    fn a_real_1ca2_contact_matches_the_golden_fixtures_bgn_end_distance_and_contact_list() {
        // VAL143 CG1 <-> VAL121 CG2, taken directly from
        // tests/fixtures/golden/1CA2.json: {"bgn": {CG1, VAL 143}, "end":
        // {CG2, VAL 121}, "distance": 4.38, "contact": ["proximal",
        // "hydrophobic"]}. Deliberately picked a HYDROPHOBIC (distance +
        // typing only) pair rather than an HBOND/WEAK_HBOND one: those
        // depend on this project's own analytically-*placed* donor
        // hydrogens (`hydrogenate.rs`), which real arpeggio instead gets
        // from OpenBabel's own placement algorithm -- confirmed directly
        // (via an ad-hoc angle check before writing this test) that a
        // nearby real GLU117 O / PHE95 CA pair is exactly this kind of
        // near-threshold case: this project's placed HA gives a WEAK_HBOND
        // angle of 133.3 degrees (passes our 130 degree cutoff), while the
        // golden fixture's real OpenBabel-derived placement apparently
        // didn't. A principled, expected divergence (see HANDOFF.md), not
        // a bug -- but the wrong case to assert byte-parity on.
        // `interacting_entities` there is "INTRA_NON_SELECTION" (that
        // fixture's real, unreproduced selection) -- not checked here, see
        // module doc.
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1CA2.cif").expect("1CA2 should load");
        let mut components = common_components();
        let zn =
            rspeggio_ccd::parser::load_ccd_component("../rspeggio-ccd/tests/fixtures/ZN_ideal.cif")
                .expect("ZN fixture should parse");
        components.insert("ZN".to_string(), zn);

        let entries = export_atom_atom_contacts(&pdb, &components)
            .expect("1CA2 is amino acids + water + zinc, all now known");

        let entry = entries
            .iter()
            .find(|e| {
                (e.bgn.label_comp_id == "VAL"
                    && e.bgn.auth_seq_id == 143
                    && e.bgn.auth_atom_id == "CG1"
                    && e.end.label_comp_id == "VAL"
                    && e.end.auth_seq_id == 121
                    && e.end.auth_atom_id == "CG2")
                    || (e.end.label_comp_id == "VAL"
                        && e.end.auth_seq_id == 143
                        && e.end.auth_atom_id == "CG1"
                        && e.bgn.label_comp_id == "VAL"
                        && e.bgn.auth_seq_id == 121
                        && e.bgn.auth_atom_id == "CG2")
            })
            .expect("VAL143 CG1 <-> VAL121 CG2 should be a real contact in 1CA2");

        assert_eq!(entry.distance, 4.38);
        assert_eq!(entry.contact, vec!["proximal", "hydrophobic"]);
    }

    #[test]
    fn a_real_1mbo_atom_identity_matches_the_golden_fixture() {
        // VAL68 CG1, taken from the atom-plane entry in
        // tests/fixtures/golden/1MBO.json used as the worked example in
        // `rings.rs`'s own design discussion -- confirms `atom_identity`'s
        // field construction (chain, seq id, comp id/type, ins code)
        // independent of the atom-atom-specific contact list.
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1MBO.cif").expect("1MBO should load");
        let mut components = common_components();
        for (comp_id, path) in [
            ("HEM", "tests/fixtures/ccd/HEM.cif"),
            ("OXY", "tests/fixtures/ccd/OXY.cif"),
        ] {
            let component = rspeggio_ccd::parser::load_ccd_component(path)
                .unwrap_or_else(|| panic!("{comp_id} fixture should parse"));
            components.insert(comp_id.to_string(), component);
        }

        let entries = export_atom_atom_contacts(&pdb, &components)
            .expect("1MBO is amino acids + water + heme + bound oxygen, all now known");

        let val68_cg1 = entries
            .iter()
            .flat_map(|e| [&e.bgn, &e.end])
            .find(|a| a.label_comp_id == "VAL" && a.auth_seq_id == 68 && a.auth_atom_id == "CG1")
            .expect("VAL68 CG1 should appear in at least one real contact");

        assert_eq!(val68_cg1.auth_asym_id, "A");
        assert_eq!(val68_cg1.label_comp_type, "P");
        assert_eq!(val68_cg1.pdbx_pdb_ins_code, " ");
    }

    #[test]
    fn a_real_water_water_contact_is_classified_correctly_in_whole_structure_mode() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");
        let components = common_components();

        let entries = export_atom_atom_contacts(&pdb, &components).expect("1UBQ should export");

        let found = entries.iter().any(|e| {
            e.bgn.label_comp_id == "HOH"
                && e.end.label_comp_id == "HOH"
                && e.interacting_entities == "WATER_WATER"
        });
        assert!(
            found,
            "expected at least one real water-water contact in 1UBQ"
        );
    }

    #[test]
    fn a_real_protein_water_contact_is_selection_water_in_whole_structure_mode() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");
        let components = common_components();

        let entries = export_atom_atom_contacts(&pdb, &components).expect("1UBQ should export");

        let found = entries.iter().any(|e| {
            (e.bgn.label_comp_id == "HOH") != (e.end.label_comp_id == "HOH")
                && e.interacting_entities == "SELECTION_WATER"
        });
        assert!(
            found,
            "expected at least one real protein-water contact in 1UBQ"
        );
    }

    #[test]
    fn a_real_protein_protein_contact_is_intra_selection_in_whole_structure_mode() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");
        let components = common_components();

        let entries = export_atom_atom_contacts(&pdb, &components).expect("1UBQ should export");

        let found = entries.iter().any(|e| {
            e.bgn.label_comp_id != "HOH"
                && e.end.label_comp_id != "HOH"
                && e.interacting_entities == "INTRA_SELECTION"
        });
        assert!(
            found,
            "expected at least one real protein-protein contact in 1UBQ"
        );
    }

    #[test]
    fn every_entry_has_exactly_one_distance_category_label() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");
        let components = common_components();

        let entries = export_atom_atom_contacts(&pdb, &components).expect("1UBQ should export");
        assert!(!entries.is_empty());

        const DISTANCE_LABELS: [&str; 5] = ["clash", "covalent", "vdw_clash", "vdw", "proximal"];
        for entry in &entries {
            let category_count = entry
                .contact
                .iter()
                .filter(|label| DISTANCE_LABELS.contains(label))
                .count();
            assert_eq!(
                category_count, 1,
                "entry {entry:?} should have exactly one distance-category label"
            );
        }
    }

    #[test]
    fn an_unknown_component_fails_the_whole_export_loudly() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");
        // Deliberately empty -- simulates every residue being unknown, the
        // same way `join.rs`'s own equivalent test does.
        let components: HashMap<String, CcdComponent> = HashMap::new();

        assert!(export_atom_atom_contacts(&pdb, &components).is_err());
    }

    #[test]
    fn the_exported_shape_actually_serializes_to_the_real_json_field_names() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");
        let components = common_components();

        let entries = export_atom_atom_contacts(&pdb, &components).expect("1UBQ should export");
        let json = serde_json::to_value(&entries[0]).expect("should serialize");

        for field in [
            "bgn",
            "end",
            "type",
            "distance",
            "contact",
            "interacting_entities",
        ] {
            assert!(
                json.get(field).is_some(),
                "expected real pdbe-arpeggio field {field:?} in serialized output"
            );
        }
        assert!(json["bgn"].get("pdbx_PDB_ins_code").is_some());
        assert_eq!(json["type"], "atom-atom");
    }
}
