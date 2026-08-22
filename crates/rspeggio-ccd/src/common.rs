// crates/rspeggio-ccd/src/common.rs
//
// The "bundled" tier of CCD lookup: a small, curated set of real CCD files
// shipped in the repo (crates/rspeggio-ccd/data/common/) covering the
// components that appear in the overwhelming majority of structures --
// standard amino acids, water, standard nucleotides, and a handful of
// common ions/crystallization additives. This resolves without any
// filesystem cache or network access.
//
// This is deliberately *not* the general-purpose CCD cache. An uncommon
// ligand or modified residue won't be here -- that's the on-disk
// cache + live-fetch + `UnknownComponent` tier, which is M4 scope, not this
// module's job.

use crate::component::CcdComponent;
use crate::parser::load_ccd_component;
use std::collections::HashMap;

const COMMON_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/data/common");

// Loads every bundled CCD file into a `comp_id -> CcdComponent` map. Uses
// `CARGO_MANIFEST_DIR` (baked in at compile time) rather than a
// cwd-relative path, so this resolves correctly regardless of where a
// consumer of this library is run from -- unlike the test fixtures
// elsewhere in this crate, which rely on `cargo test`'s working directory.
//
// A missing bundle directory or a bundled file that fails to parse is a
// packaging defect in this crate, not a runtime condition a caller can
// meaningfully handle, so both panic rather than returning a partial map.
pub fn common_components() -> HashMap<String, CcdComponent> {
    let entries = std::fs::read_dir(COMMON_DIR)
        .unwrap_or_else(|e| panic!("bundled common CCD directory missing: {COMMON_DIR}: {e}"));

    let mut components = HashMap::new();
    for entry in entries {
        let path = entry
            .unwrap_or_else(|e| panic!("failed to read bundled CCD directory entry: {e}"))
            .path();
        if path.extension().and_then(|e| e.to_str()) != Some("cif") {
            continue;
        }
        let comp_id = path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or_else(|| panic!("bundled CCD file has no usable name: {path:?}"))
            .to_string();
        let path_str = path
            .to_str()
            .unwrap_or_else(|| panic!("bundled CCD file path is not valid UTF-8: {path:?}"));
        let component = load_ccd_component(path_str)
            .unwrap_or_else(|| panic!("bundled common CCD file failed to parse: {path_str}"));
        components.insert(comp_id, component);
    }
    components
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bundles_all_expected_common_components() {
        let components = common_components();

        assert_eq!(components.len(), 37);

        for amino_acid in [
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS",
            "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        ] {
            assert!(
                components.contains_key(amino_acid),
                "missing standard amino acid {amino_acid}"
            );
        }

        for other in [
            "HOH", "DA", "DC", "DG", "DT", "A", "C", "G", "U", "NA", "MG", "CA", "K", "CL", "SO4",
            "PO4", "GOL",
        ] {
            assert!(
                components.contains_key(other),
                "missing common component {other}"
            );
        }
    }

    #[test]
    fn bundled_alanine_has_the_expected_shape() {
        let components = common_components();
        let ala = components.get("ALA").expect("ALA should be bundled");

        assert_eq!(
            ala.atoms().len(),
            13,
            "ALA: N,CA,C,O,CB + 8 H (ideal coords)"
        );
        assert!(!ala.bonds().is_empty());
    }

    #[test]
    fn bundled_water_has_two_hydrogens_bonded_to_the_oxygen() {
        let components = common_components();
        let hoh = components.get("HOH").expect("HOH should be bundled");

        assert_eq!(hoh.atoms().len(), 3, "O, H1, H2");
        assert_eq!(hoh.bonds().len(), 2, "O-H1, O-H2");
    }

    #[test]
    fn bundled_calcium_ion_is_a_single_bondless_atom() {
        let components = common_components();
        // CA's own CCD file uses CIF's scalar (non-loop) single-row form,
        // same as ZN in rspeggio_ccd's parser tests -- this exercises that
        // path against a real bundled file, not just a fixture.
        let ca = components.get("CA").expect("CA should be bundled");

        assert_eq!(ca.atoms().len(), 1);
        assert!(ca.bonds().is_empty());
    }
}
