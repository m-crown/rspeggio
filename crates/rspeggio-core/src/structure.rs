// crates/rspeggio-core/src/structure.rs

#[derive(Debug, PartialEq)]
pub struct StructureAtom {
    pub serial: usize,
    pub name: String,
    pub element: String,
    pub pos: (f64, f64, f64),
}

#[derive(Debug, PartialEq)]
pub struct StructureResidue {
    pub comp_id: String,
    pub chain_id: String,
    pub seq_id: isize,
    pub atoms: Vec<StructureAtom>,
}

pub struct Structure {
    pub residues: Vec<StructureResidue>,
}

impl Structure {
    // Flattens a parsed `pdbtbx::PDB` into rspeggio's own minimal
    // residue/atom shape. Doesn't touch CCD data yet -- this is just the
    // structural geometry, joining each residue against its `CcdComponent`
    // is the next step once this exists.
    pub fn from_pdb(pdb: &pdbtbx::PDB) -> Self {
        let residues = pdb
            .chains()
            .flat_map(|chain| {
                let chain_id = chain.id().to_string();
                chain.residues().map(move |residue| StructureResidue {
                    comp_id: residue.name().unwrap_or_default().to_string(),
                    chain_id: chain_id.clone(),
                    seq_id: residue.id().0,
                    atoms: residue
                        .atoms()
                        .map(|atom| StructureAtom {
                            serial: atom.serial_number(),
                            name: atom.name().to_string(),
                            element: atom
                                .element()
                                .map(|e| e.symbol().to_string())
                                .unwrap_or_default(),
                            pos: atom.pos(),
                        })
                        .collect(),
                })
            })
            .collect();
        Structure { residues }
    }
}

// Opens and converts a structure file from disk in one step. `None` covers
// both an unreadable/unparseable file and a fatal parse error at
// `StrictnessLevel::Medium` -- callers needing to distinguish those should
// call `pdbtbx::open`/`Structure::from_pdb` directly instead.
pub fn load_structure(path: &str) -> Option<Structure> {
    let (pdb, errors) = pdbtbx::open(path).ok()?;
    if errors
        .iter()
        .any(|e| e.fails(pdbtbx::StrictnessLevel::Medium))
    {
        return None;
    }
    Some(Structure::from_pdb(&pdb))
}

#[cfg(test)]
mod tests {
    use super::*;
    use pdbtbx::StrictnessLevel;

    #[test]
    fn opens_a_real_mmcif_structure() {
        let (pdb, errors) = pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif")
            .expect("1UBQ fixture should parse");

        assert!(
            errors.iter().all(|e| !e.fails(StrictnessLevel::Medium)),
            "no fatal parse errors: {errors:?}"
        );
        assert!(pdb.atom_count() > 0);
        assert!(pdb.residue_count() > 0);
    }

    #[test]
    fn reads_the_first_residue_of_ubiquitin() {
        let (pdb, _errors) = pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif")
            .expect("1UBQ fixture should parse");

        let first_residue = pdb
            .residues()
            .find(|r| r.name() == Some("MET"))
            .expect("ubiquitin starts with MET1");

        assert_eq!(first_residue.serial_number(), 1);
        assert!(first_residue.atoms().any(|a| a.name() == "CA"));
    }

    #[test]
    fn ubiquitin_fixture_includes_waters_as_a_separate_hetero_residue() {
        let (pdb, _errors) = pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif")
            .expect("1UBQ fixture should parse");

        assert!(
            pdb.residues().any(|r| r.name() == Some("HOH")),
            "fixture is protein + waters per the plan's fixture table"
        );
    }

    #[test]
    fn converts_a_real_pdb_into_the_internal_structure_shape() {
        let (pdb, _errors) = pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif")
            .expect("1UBQ fixture should parse");

        let structure = Structure::from_pdb(&pdb);

        assert_eq!(structure.residues.len(), pdb.residue_count());

        let met1 = structure
            .residues
            .iter()
            .find(|r| r.comp_id == "MET" && r.seq_id == 1)
            .expect("MET1 should be present");

        assert_eq!(met1.chain_id, "A");
        let ca = met1
            .atoms
            .iter()
            .find(|a| a.name == "CA")
            .expect("MET1 should have a CA atom");
        assert_eq!(ca.element, "C");
        assert_ne!(ca.pos, (0.0, 0.0, 0.0));
    }

    #[test]
    fn loads_a_structure_straight_from_a_file_path() {
        let structure =
            load_structure("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");

        assert!(!structure.residues.is_empty());
    }

    #[test]
    fn load_structure_returns_none_for_a_missing_file() {
        assert!(load_structure("../../tests/fixtures/structures/does_not_exist.cif").is_none());
    }
}
