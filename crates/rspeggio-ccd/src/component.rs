#[derive(Debug, PartialEq)]
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

#[derive(Debug, PartialEq)]
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

#[derive(Debug, PartialEq)]
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

pub struct CcdComponent {
    atoms: Vec<CcdAtom>,
    bonds: Vec<CcdBond>,
}

impl CcdComponent {
    pub fn new(atoms: Vec<CcdAtom>, bonds: Vec<CcdBond>) -> Self {
        Self { atoms, bonds }
    }

    pub fn atoms(&self) -> &[CcdAtom] {
        &self.atoms
    }

    pub fn bonds(&self) -> &[CcdBond] {
        &self.bonds
    }
}
