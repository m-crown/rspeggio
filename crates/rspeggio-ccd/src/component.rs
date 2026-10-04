#[derive(Debug, Clone, PartialEq)]
pub struct CcdAtom {
    pub(crate) atom_id: String,
    pub(crate) element: String,
    pub(crate) aromatic: bool,
    // From `_chem_comp_atom.pdbx_leaving_atom_flag`: true for an atom that
    // exists only in this component's free/monomeric form and is removed
    // when it polymerizes into a chain in a non-terminal position (e.g. a
    // standard amino acid's OXT/HXT and the second amine hydrogen). The
    // atom that survives polymerization (e.g. the backbone amide H itself)
    // is never a leaving atom -- that invariant is what lets binary
    // donor/acceptor typing ignore this distinction; only a count of
    // donor hydrogens needs to exclude leaving atoms for internal residues.
    pub(crate) leaving: bool,
}

impl CcdAtom {
    pub fn new(atom_id: String, element: String, aromatic: bool, leaving: bool) -> Self {
        Self {
            atom_id,
            element,
            aromatic,
            leaving,
        }
    }

    pub fn atom_id(&self) -> &str {
        &self.atom_id
    }

    pub fn element(&self) -> &str {
        &self.element
    }

    pub fn aromatic(&self) -> bool {
        self.aromatic
    }

    pub fn leaving(&self) -> bool {
        self.leaving
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Aromatic,
}

impl BondOrder {
    // Maps the raw `value_order` token (`SING`/`DOUB`/`TRIP`/`AROM`) from a
    // `_chem_comp_bond` row. `None` for anything else, so an unrecognized
    // value fails the row rather than silently defaulting to some order.
    pub fn from_ccd_str(s: &str) -> Option<Self> {
        match s {
            "SING" => Some(Self::Single),
            "DOUB" => Some(Self::Double),
            "TRIP" => Some(Self::Triple),
            "AROM" => Some(Self::Aromatic),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct CcdBond {
    pub(crate) atom_id_1: String,
    pub(crate) atom_id_2: String,
    pub(crate) order: BondOrder,
    pub(crate) aromatic: bool,
}

impl CcdBond {
    pub fn new(atom_id_1: String, atom_id_2: String, order: BondOrder, aromatic: bool) -> Self {
        Self {
            atom_id_1,
            atom_id_2,
            order,
            aromatic,
        }
    }

    pub fn atom_id_1(&self) -> &str {
        &self.atom_id_1
    }

    pub fn atom_id_2(&self) -> &str {
        &self.atom_id_2
    }

    pub fn order(&self) -> &BondOrder {
        &self.order
    }

    pub fn aromatic(&self) -> bool {
        self.aromatic
    }
}

// Real pdbe-arpeggio's `ComponentType` (`config.py`): the classification
// exported as JSON's `label_comp_type`, one letter per component. Derived
// from the CCD's own `_chem_comp.type` (a real, small controlled
// vocabulary -- polymer-linking type by polymer kind, or "NON-POLYMER"/
// "OTHER"), with one further real refinement: a `NON-POLYMER` component
// splits into `Water` or `BoundMolecule` by checking whether its
// `_chem_comp.name` is literally "WATER" (there's no dedicated water type
// in the CCD vocabulary itself).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ComponentType {
    Polypeptide,
    PolyDeoxyribonucleotide,
    PolyRibonucleotide,
    Polysaccharide,
    Other,
    BoundMolecule,
    Water,
}

impl ComponentType {
    // The one-letter code real pdbe-arpeggio's JSON export actually uses
    // for `label_comp_type` (`config.ComponentType`'s enum member names).
    pub fn code(&self) -> &'static str {
        match self {
            Self::Polypeptide => "P",
            Self::PolyDeoxyribonucleotide => "D",
            Self::PolyRibonucleotide => "R",
            Self::Polysaccharide => "S",
            Self::Other => "O",
            Self::BoundMolecule => "B",
            Self::Water => "W",
        }
    }

    // Ported verbatim from real pdbe-arpeggio's `ComponentType.from_chem_comp_type`
    // lookup table (`config.py`) -- every raw `_chem_comp.type` string the
    // real CCD vocabulary actually uses, matched case-insensitively (the
    // real code uppercases first). `None` for anything not in that table,
    // same failure posture as `BondOrder::from_ccd_str`: an unrecognized
    // value fails the parse rather than guessing a type.
    pub fn from_chem_comp_type(raw_type: &str, name: &str) -> Option<Self> {
        match raw_type.to_uppercase().as_str() {
            "D-BETA-PEPTIDE, C-GAMMA LINKING"
            | "D-GAMMA-PEPTIDE, C-DELTA LINKING"
            | "D-PEPTIDE COOH CARBOXY TERMINUS"
            | "D-PEPTIDE NH3 AMINO TERMINUS"
            | "D-PEPTIDE LINKING"
            | "L-BETA-PEPTIDE, C-GAMMA LINKING"
            | "L-GAMMA-PEPTIDE, C-DELTA LINKING"
            | "L-PEPTIDE COOH CARBOXY TERMINUS"
            | "L-PEPTIDE NH3 AMINO TERMINUS"
            | "L-PEPTIDE LINKING"
            | "PEPTIDE LINKING"
            | "PEPTIDE-LIKE" => Some(Self::Polypeptide),
            "DNA OH 3 PRIME TERMINUS"
            | "DNA OH 5 PRIME TERMINUS"
            | "DNA LINKING"
            | "L-DNA LINKING" => Some(Self::PolyDeoxyribonucleotide),
            "L-RNA LINKING"
            | "RNA OH 3 PRIME TERMINUS"
            | "RNA OH 5 PRIME TERMINUS"
            | "RNA LINKING" => Some(Self::PolyRibonucleotide),
            "D-SACCHARIDE"
            | "D-SACCHARIDE, ALPHA LINKING"
            | "D-SACCHARIDE, BETA LINKING"
            | "L-SACCHARIDE"
            | "L-SACCHARIDE, ALPHA LINKING"
            | "L-SACCHARIDE, BETA LINKING"
            | "SACCHARIDE" => Some(Self::Polysaccharide),
            "OTHER" => Some(Self::Other),
            "NON-POLYMER" if name.eq_ignore_ascii_case("WATER") => Some(Self::Water),
            "NON-POLYMER" => Some(Self::BoundMolecule),
            _ => None,
        }
    }
}

#[derive(Debug, Clone)]
pub struct CcdComponent {
    atoms: Vec<CcdAtom>,
    bonds: Vec<CcdBond>,
    component_type: ComponentType,
}

impl CcdComponent {
    pub fn new(atoms: Vec<CcdAtom>, bonds: Vec<CcdBond>, component_type: ComponentType) -> Self {
        Self {
            atoms,
            bonds,
            component_type,
        }
    }

    pub fn atoms(&self) -> &[CcdAtom] {
        &self.atoms
    }

    pub fn bonds(&self) -> &[CcdBond] {
        &self.bonds
    }

    pub fn component_type(&self) -> ComponentType {
        self.component_type
    }
}
