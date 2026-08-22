// crates/rspeggio-core/src/structure.rs

#[cfg(test)]
mod tests {
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
}
