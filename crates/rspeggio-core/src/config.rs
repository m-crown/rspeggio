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
// Source: Open Babel's real `elementtable.h`
// (https://github.com/openbabel/openbabel/blob/master/src/elementtable.h),
// itself compiled from the Blue Obelisk Cheminformatics Data Repository
// (http://www.blueobelisk.org/repos/blueobelisk/elements.xml). Real
// pdbe-arpeggio does not hardcode these values either -- it pulls the
// same numbers at runtime via `ob.GetCovalentRad`/`ob.GetVdwRad`
// (interactions.py:1501/1509 in PDBeurope/arpeggio). Porting the published
// numbers directly, rather than linking Open Babel, keeps decision 02's
// no-FFI rule intact while still matching the oracle's actual values.
//
// The full real table (all 118 elements H through Og, fetched directly
// from the source above, not hand-picked) -- ported completely rather than
// scoped to the current fixture corpus, so a future structure/ligand with
// an element this project hasn't seen yet doesn't silently fall through to
// `None`. A handful of the heaviest synthetic elements share Open Babel's
// own literal placeholder values (1.60/2.00 -- documented in its own file
// header as "if unknown"), not real measurements; ported verbatim anyway,
// since that's genuinely what the oracle itself uses for them too.
const COVALENT_RADII: &[(&str, f64)] = &[
    ("H", 0.31),
    ("HE", 0.28),
    ("LI", 1.28),
    ("BE", 0.96),
    ("B", 0.84),
    ("C", 0.76),
    ("N", 0.71),
    ("O", 0.66),
    ("F", 0.57),
    ("NE", 0.58),
    ("NA", 1.66),
    ("MG", 1.41),
    ("AL", 1.21),
    ("SI", 1.11),
    ("P", 1.07),
    ("S", 1.05),
    ("CL", 1.02),
    ("AR", 1.06),
    ("K", 2.03),
    ("CA", 1.76),
    ("SC", 1.7),
    ("TI", 1.6),
    ("V", 1.53),
    ("CR", 1.39),
    ("MN", 1.39),
    ("FE", 1.32),
    ("CO", 1.26),
    ("NI", 1.24),
    ("CU", 1.32),
    ("ZN", 1.22),
    ("GA", 1.22),
    ("GE", 1.2),
    ("AS", 1.19),
    ("SE", 1.2),
    ("BR", 1.2),
    ("KR", 1.16),
    ("RB", 2.2),
    ("SR", 1.95),
    ("Y", 1.9),
    ("ZR", 1.75),
    ("NB", 1.64),
    ("MO", 1.54),
    ("TC", 1.47),
    ("RU", 1.46),
    ("RH", 1.42),
    ("PD", 1.39),
    ("AG", 1.45),
    ("CD", 1.44),
    ("IN", 1.42),
    ("SN", 1.39),
    ("SB", 1.39),
    ("TE", 1.38),
    ("I", 1.39),
    ("XE", 1.4),
    ("CS", 2.44),
    ("BA", 2.15),
    ("LA", 2.07),
    ("CE", 2.04),
    ("PR", 2.03),
    ("ND", 2.01),
    ("PM", 1.99),
    ("SM", 1.98),
    ("EU", 1.98),
    ("GD", 1.96),
    ("TB", 1.94),
    ("DY", 1.92),
    ("HO", 1.92),
    ("ER", 1.89),
    ("TM", 1.9),
    ("YB", 1.87),
    ("LU", 1.87),
    ("HF", 1.75),
    ("TA", 1.7),
    ("W", 1.62),
    ("RE", 1.51),
    ("OS", 1.44),
    ("IR", 1.41),
    ("PT", 1.36),
    ("AU", 1.36),
    ("HG", 1.32),
    ("TL", 1.45),
    ("PB", 1.46),
    ("BI", 1.48),
    ("PO", 1.4),
    ("AT", 1.5),
    ("RN", 1.5),
    ("FR", 2.6),
    ("RA", 2.21),
    ("AC", 2.15),
    ("TH", 2.06),
    ("PA", 2.0),
    ("U", 1.96),
    ("NP", 1.9),
    ("PU", 1.87),
    ("AM", 1.8),
    ("CM", 1.69),
    ("BK", 1.6),
    ("CF", 1.6),
    ("ES", 1.6),
    ("FM", 1.6),
    ("MD", 1.6),
    ("NO", 1.6),
    ("LR", 1.6),
    ("RF", 1.6),
    ("DB", 1.6),
    ("SG", 1.6),
    ("BH", 1.6),
    ("HS", 1.6),
    ("MT", 1.6),
    ("DS", 1.6),
    ("RG", 1.6),
    ("CN", 1.6),
    ("NH", 1.6),
    ("FL", 1.6),
    ("MC", 1.6),
    ("LV", 1.6),
    ("TS", 1.6),
    ("OG", 1.6),
];

const VDW_RADII: &[(&str, f64)] = &[
    ("H", 1.1),
    ("HE", 1.4),
    ("LI", 1.81),
    ("BE", 1.53),
    ("B", 1.92),
    ("C", 1.7),
    ("N", 1.55),
    ("O", 1.52),
    ("F", 1.47),
    ("NE", 1.54),
    ("NA", 2.27),
    ("MG", 1.73),
    ("AL", 1.84),
    ("SI", 2.1),
    ("P", 1.8),
    ("S", 1.8),
    ("CL", 1.75),
    ("AR", 1.88),
    ("K", 2.75),
    ("CA", 2.31),
    ("SC", 2.3),
    ("TI", 2.15),
    ("V", 2.05),
    ("CR", 2.05),
    ("MN", 2.05),
    ("FE", 2.05),
    ("CO", 2.0),
    ("NI", 2.0),
    ("CU", 2.0),
    ("ZN", 2.1),
    ("GA", 1.87),
    ("GE", 2.11),
    ("AS", 1.85),
    ("SE", 1.9),
    ("BR", 1.83),
    ("KR", 2.02),
    ("RB", 3.03),
    ("SR", 2.49),
    ("Y", 2.4),
    ("ZR", 2.3),
    ("NB", 2.15),
    ("MO", 2.1),
    ("TC", 2.05),
    ("RU", 2.05),
    ("RH", 2.0),
    ("PD", 2.05),
    ("AG", 2.1),
    ("CD", 2.2),
    ("IN", 2.2),
    ("SN", 1.93),
    ("SB", 2.17),
    ("TE", 2.06),
    ("I", 1.98),
    ("XE", 2.16),
    ("CS", 3.43),
    ("BA", 2.68),
    ("LA", 2.5),
    ("CE", 2.48),
    ("PR", 2.47),
    ("ND", 2.45),
    ("PM", 2.43),
    ("SM", 2.42),
    ("EU", 2.4),
    ("GD", 2.38),
    ("TB", 2.37),
    ("DY", 2.35),
    ("HO", 2.33),
    ("ER", 2.32),
    ("TM", 2.3),
    ("YB", 2.28),
    ("LU", 2.27),
    ("HF", 2.25),
    ("TA", 2.2),
    ("W", 2.1),
    ("RE", 2.05),
    ("OS", 2.0),
    ("IR", 2.0),
    ("PT", 2.05),
    ("AU", 2.1),
    ("HG", 2.05),
    ("TL", 1.96),
    ("PB", 2.02),
    ("BI", 2.07),
    ("PO", 1.97),
    ("AT", 2.02),
    ("RN", 2.2),
    ("FR", 3.48),
    ("RA", 2.83),
    ("AC", 2.0),
    ("TH", 2.4),
    ("PA", 2.0),
    ("U", 2.3),
    ("NP", 2.0),
    ("PU", 2.0),
    ("AM", 2.0),
    ("CM", 2.0),
    ("BK", 2.0),
    ("CF", 2.0),
    ("ES", 2.0),
    ("FM", 2.0),
    ("MD", 2.0),
    ("NO", 2.0),
    ("LR", 2.0),
    ("RF", 2.0),
    ("DB", 2.0),
    ("SG", 2.0),
    ("BH", 2.0),
    ("HS", 2.0),
    ("MT", 2.0),
    ("DS", 2.0),
    ("RG", 2.0),
    ("CN", 2.0),
    ("NH", 2.0),
    ("FL", 2.0),
    ("MC", 2.0),
    ("LV", 2.0),
    ("TS", 2.0),
    ("OG", 2.0),
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
    fn the_full_periodic_table_is_now_covered_not_just_the_original_fixture_corpus() {
        // Selenium (real: selenomethionine SAD-phasing substitutes this
        // for sulfur in real crystallography, so a real structure could
        // plausibly contain it) and Gold (a real heavy-atom derivative
        // element) are both genuine elements from real chemistry that the
        // table's original, narrower scope (H/C/N/O/F/Na/Mg/P/S/Cl/K/Ca/
        // Mn/Fe/Co/Ni/Cu/Zn/Br/I) didn't cover -- confirming they resolve
        // now is the actual regression test for the scope expansion, not
        // just re-checking an element already covered before.
        assert_eq!(covalent_radius("SE"), Some(1.20));
        assert_eq!(vdw_radius("SE"), Some(1.90));
        assert_eq!(covalent_radius("AU"), Some(1.36));
        assert_eq!(vdw_radius("AU"), Some(2.10));
    }

    #[test]
    fn every_real_element_has_both_a_covalent_and_vdw_radius() {
        assert_eq!(
            COVALENT_RADII.len(),
            118,
            "the real periodic table (H through Og) has 118 elements"
        );
        assert_eq!(VDW_RADII.len(), 118);
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
