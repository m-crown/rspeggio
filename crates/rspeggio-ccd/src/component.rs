// temporary remove once used downstream
#![allow(dead_code)]

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
}

pub struct CcdComponent {
    atoms: Vec<CcdAtom>,
}
