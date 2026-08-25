// enabled whilst structs are not in use. remove before finalisation
#![allow(dead_code)]

use bitflags::bitflags;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DistanceCategory {
    Clash,
    Covalent,
    VdwClash,
    Vdw,
    Proximal,
} // positions 0-4 in SIFt representation see pdbe-arpeggio `interactions.py:748-773`

// Compensation factor added to the VdW-radii-sum upper bound before an
// atom pair falls through to Proximal. Matches real pdbe-arpeggio's
// default `vdw_comp` parameter (`interactions.py:37`).
pub const VDW_COMP_FACTOR: f64 = 0.1;

// u16 bit shift to build the FeatureBits e.g. HBOND = 0000000000000001 + WEAK_HBOND = 0000000000000010
bitflags! {
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub struct FeatureBits: u16 {
        const HBOND         = 1 << 0;
        const WEAK_HBOND    = 1 << 1;
        const XBOND         = 1 << 2;
        const IONIC         = 1 << 3;
        const METAL_COMPLEX = 1 << 4;
        const AROMATIC      = 1 << 5;
        const HYDROPHOBIC   = 1 << 6;
        const CARBONYL      = 1 << 7;
        const POLAR         = 1 << 8;
        const WEAK_POLAR    = 1 << 9;
    }
}

struct Sift {
    distance_category: DistanceCategory,
    features: FeatureBits,
}

pub struct HbondThresholds {
    pub distance: f64,
    pub polar_distance: f64,
    pub angle_degrees: f64,
}

pub struct WeakHbondThresholds {
    pub distance: f64,
    pub weak_polar_distance: f64,
    pub angle_degrees: f64,
}

pub struct IonicThresholds {
    pub distance: f64,
}

pub struct HydrophobicThresholds {
    pub distance: f64,
}

pub struct CarbonylThresholds {
    pub distance: f64,
}

pub struct MetalThresholds {
    pub distance: f64,
}

pub struct AromaticThresholds {
    pub distance: f64,
}

pub struct XbondThresholds {
    // Real pdbe-arpeggio's `angle theta 1` (`interactions.py`'s `is_xbond`):
    // the C-X...acceptor angle must be at least this wide (a loose,
    // one-sided lower bound only -- real arpeggio's config also defines a
    // `theta 2 min/max` range and a `catmap distance`, but neither is ever
    // actually read anywhere in its own codebase, confirmed by searching
    // it, so this only ports the threshold that's genuinely load-bearing).
    pub angle_theta_1_degrees: f64,
}

pub const CONTACT_TYPES_MAX_DIST: f64 = 4.5;

pub const HBOND: HbondThresholds = HbondThresholds {
    distance: 3.9,
    polar_distance: 3.5,
    angle_degrees: 90.0,
};

pub const WEAK_HBOND: WeakHbondThresholds = WeakHbondThresholds {
    distance: 3.6,
    weak_polar_distance: 3.5,
    angle_degrees: 130.0,
};

pub const IONIC: IonicThresholds = IonicThresholds { distance: 4.0 };

pub const HYDROPHOBIC: HydrophobicThresholds = HydrophobicThresholds { distance: 4.5 };

pub const CARBONYL: CarbonylThresholds = CarbonylThresholds { distance: 3.6 };

pub const METAL: MetalThresholds = MetalThresholds { distance: 2.8 };

// Simple atom-atom aromatic contact: both atoms typed aromatic and within
// this distance (real pdbe-arpeggio's `CONTACT_TYPES['aromatic']['distance']`,
// `interactions.py:916`) -- distinct from `rings::CENTROID_DISTANCE_MAX`,
// which gates the separate ring-*plane*-geometry (centroid/normal)
// contacts in `rings.rs`, not this per-atom-pair feature bit.
pub const AROMATIC: AromaticThresholds = AromaticThresholds { distance: 4.0 };

pub const XBOND: XbondThresholds = XbondThresholds {
    angle_theta_1_degrees: 120.0,
};

pub const MAINCHAIN_ATOMS: [&str; 5] = ["N", "C", "CA", "O", "OXT"];

// The 20 standard amino acids, ported verbatim from real pdbe-arpeggio's
// `STD_RES` (config.py). Used to decide which residues participate in
// polypeptide-chain sequence-adjacency (waters, ions, and ligands never
// do, regardless of numbering).
pub const STANDARD_AMINO_ACIDS: [&str; 20] = [
    "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO",
    "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR",
];

// Covalent and van der Waals radii (Å), keyed by element symbol (matched
// case-insensitively against the CCD's `type_symbol`, which is uppercase).
//
// Source: Open Babel's `elementtable.h`
// (https://github.com/openbabel/openbabel/blob/master/src/elementtable.h),
// itself compiled from the Blue Obelisk Cheminformatics Data Repository
// (http://www.blueobelisk.org/repos/blueobelisk/elements.xml). Real
// pdbe-arpeggio does not hardcode these values either -- it pulls the
// same numbers at runtime via `ob.GetCovalentRad`/`ob.GetVdwRad`
// (interactions.py:1501/1509 in PDBeurope/arpeggio). Porting the published
// numbers directly, rather than linking Open Babel, keeps decision 02's
// no-FFI rule intact while still matching the oracle's actual values.
//
// TODO: revisit and extend as real fixture structures introduce elements
// not covered here. Currently limited to what the corpus + common bundle
// actually reference: H, C, N, O, F, Na, Mg, P, S, Cl, K, Ca, Mn, Fe, Co,
// Ni, Cu, Zn, Br, I.
const COVALENT_RADII: &[(&str, f64)] = &[
    ("H", 0.31),
    ("C", 0.76),
    ("N", 0.71),
    ("O", 0.66),
    ("F", 0.57),
    ("NA", 1.66),
    ("MG", 1.41),
    ("P", 1.07),
    ("S", 1.05),
    ("CL", 1.02),
    ("K", 2.03),
    ("CA", 1.76),
    ("MN", 1.39),
    ("FE", 1.32),
    ("CO", 1.26),
    ("NI", 1.24),
    ("CU", 1.32),
    ("ZN", 1.22),
    ("BR", 1.20),
    ("I", 1.39),
];

const VDW_RADII: &[(&str, f64)] = &[
    ("H", 1.10),
    ("C", 1.70),
    ("N", 1.55),
    ("O", 1.52),
    ("F", 1.47),
    ("NA", 2.27),
    ("MG", 1.73),
    ("P", 1.80),
    ("S", 1.80),
    ("CL", 1.75),
    ("K", 2.75),
    ("CA", 2.31),
    ("MN", 2.05),
    ("FE", 2.05),
    ("CO", 2.00),
    ("NI", 2.00),
    ("CU", 2.00),
    ("ZN", 2.10),
    ("BR", 1.83),
    ("I", 1.98),
];

pub fn covalent_radius(element: &str) -> Option<f64> {
    COVALENT_RADII
        .iter()
        .find(|(symbol, _)| symbol.eq_ignore_ascii_case(element))
        .map(|(_, radius)| *radius)
}

pub fn vdw_radius(element: &str) -> Option<f64> {
    VDW_RADII
        .iter()
        .find(|(symbol, _)| symbol.eq_ignore_ascii_case(element))
        .map(|(_, radius)| *radius)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn combining_flags_works() {
        let f = FeatureBits::HBOND | FeatureBits::AROMATIC;
        assert!(f.contains(FeatureBits::HBOND));
        assert!(f.contains(FeatureBits::AROMATIC));
        assert!(!f.contains(FeatureBits::IONIC));
    }

    #[test]
    fn validate_thresholds() {
        assert_eq!(HBOND.distance, 3.9);
        assert_eq!(IONIC.distance, 4.0);
        assert_eq!(HYDROPHOBIC.distance, 4.5);
        assert_eq!(CARBONYL.distance, 3.6);
        assert_eq!(METAL.distance, 2.8);
        assert_eq!(CONTACT_TYPES_MAX_DIST, 4.5);
    }

    #[test]
    fn validate_atoms() {
        assert!(MAINCHAIN_ATOMS.contains(&"CA"));
        assert!(!MAINCHAIN_ATOMS.contains(&"CB"));
    }

    #[test]
    fn radii_match_openbabels_published_values() {
        // Spot-check against the real elementtable.h values, not just
        // whatever this table happens to already say.
        assert_eq!(covalent_radius("C"), Some(0.76));
        assert_eq!(vdw_radius("C"), Some(1.70));
        assert_eq!(covalent_radius("O"), Some(0.66));
        assert_eq!(vdw_radius("O"), Some(1.52));
        assert_eq!(covalent_radius("ZN"), Some(1.22));
        assert_eq!(vdw_radius("ZN"), Some(2.10));
    }

    #[test]
    fn radii_lookup_is_case_insensitive() {
        // CCD type_symbol is uppercase ("ZN"), but callers shouldn't have
        // to know or care about that convention.
        assert_eq!(covalent_radius("zn"), covalent_radius("ZN"));
        assert_eq!(vdw_radius("Cl"), vdw_radius("CL"));
    }

    #[test]
    fn radii_are_none_for_an_unlisted_element() {
        assert_eq!(covalent_radius("XX"), None);
        assert_eq!(vdw_radius("XX"), None);
    }

    #[test]
    fn every_covalent_entry_has_a_matching_vdw_entry() {
        // The two tables are meant to describe the same element set --
        // catches a copy-paste gap between them if one is extended and
        // the other forgotten.
        for (symbol, _) in COVALENT_RADII {
            assert!(
                vdw_radius(symbol).is_some(),
                "{symbol} has a covalent radius but no vdw radius"
            );
        }
        for (symbol, _) in VDW_RADII {
            assert!(
                covalent_radius(symbol).is_some(),
                "{symbol} has a vdw radius but no covalent radius"
            );
        }
    }
}
