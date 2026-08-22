#[derive(Debug, PartialEq)]
pub struct CcdAtom {
    pub(crate) atom_id: String,
    pub(crate) element: String,
    pub(crate) aromatic: bool,
}

impl CcdAtom {
    pub fn new(atom_id: String, element: String, aromatic: bool) -> Self {
        Self {
            atom_id,
            element,
            aromatic,
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
