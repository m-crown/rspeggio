// temporary remove once used downstream
#![allow(dead_code)]

pub struct CcdAtom {
    atom_id: String,
    element: String,
    aromatic: bool,
}

pub struct CcdComponent {
    atoms: Vec<CcdAtom>,
}
