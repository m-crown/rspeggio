// enabled whilst structs are not in use. remove before finalisation
#![allow(dead_code)]

use bitflags::bitflags;

enum DistanceCategory {
    Clash,
    Covalent,
    VdwClash,
    Vdw,
    Proximal,
} // positions 0-4 in SIFt representation see pdbe-arpeggio `interactions.py:748-773`

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
    distance: f64,
    polar_distance: f64,
    angle_degrees: f64,
}

pub struct WeakHbondThresholds {
    distance: f64,
    weak_polar_distance: f64,
    angle_degrees: f64,
}

pub struct IonicThresholds {
    distance: f64,
}

pub struct HydrophobicThresholds {
    distance: f64,
}

pub struct CarbonylThresholds {
    distance: f64,
}

pub struct MetalThresholds {
    distance: f64,
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

pub const MAINCHAIN_ATOMS: [&str; 5] = ["N", "C", "CA", "O", "OXT"];

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
}
