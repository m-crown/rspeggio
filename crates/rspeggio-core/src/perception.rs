// crates/rspeggio-core/src/perception.rs
//
// M10: geometric-perception fallback for components with no CCD entry.
// The default pipeline is CCD-first and fails loudly on an unknown
// comp_id (decision 03, `export::component_for`). This module is the
// opt-in alternative: perceive a component's bond graph straight from one
// real residue instance's 3D coordinates, producing an ordinary
// `CcdComponent` that every downstream consumer (typing, features, rings,
// export) already accepts unchanged.
//
// Nothing here runs unless a caller asks for it:
// `augment_with_geometric_fallback` fills the gaps in a components map
// before it's passed to the export functions, whose signatures are
// untouched. Without that call, unknown components still fail loudly.
//
// Four passes, each checked against real CCD entries in the tests below:
//   1. Connectivity: two atoms in the same residue are bonded if
//      0.4 A < d < summed covalent radii + 0.45 A -- Open Babel's own
//      `ConnectTheDots` window, the rule real pdbe-arpeggio's bond graph
//      comes from. Bonds never cross a residue boundary, matching what a
//      real CCD entry contains.
//   2. Initial bond order: `Double` if d < `config::BOND_ORDER_DOUBLE_RATIO`
//      (0.94) x the summed radii, else `Single` (calibrated against real
//      1UBQ bond lengths, see that constant's doc). Resonance-shortened
//      single bonds -- amide C-N, both carboxylate C-O, guanidinium C-N --
//      land on the `Double` side here; distance alone can't tell them
//      apart from real double bonds.
//   3. Valence repair: while any C/N/O has a bond-order sum above its
//      standard valence (4/3/2), demote its longest `Double` bond to
//      `Single`. This recovers the CCD's own single-Kekule-form drawing of
//      those resonance groups: a carboxylate carbon (R + two "double" C-O =
//      5) keeps the shorter C=O; an amide carbon (R + C=O + "double" C-N =
//      5) demotes the longer C-N; a guanidinium carbon (three "double" C-N
//      = 6) keeps one. Same idea as Open Babel's own valence-driven
//      cleanup, reduced to the one rule these groups need.
//   4. Implicit hydrogens: most X-ray structures deposit no hydrogens, but
//      typing decides donor status from H neighbours in the component's
//      graph (a CCD entry lists every H whether or not it was deposited).
//      Each C/N/O/S gets (standard valence - bond-order sum) hydrogens
//      (C 4, N 3, O 2, S 2), added as atoms with no coordinates and named
//      `<parent>_H<n>` so they can't collide with a real atom name.
//      `features.rs` then places them analytically exactly as it does for
//      a CCD hydrogen the structure didn't resolve. This is the role Open
//      Babel's `AddHydrogens` plays in real pdbe-arpeggio.
//   Component type is always `BoundMolecule` (the real CCD classification
//   for a non-water `NON-POLYMER`); polymer type isn't recoverable from
//   one residue's shape.
//
// Known limitations (deliberate v1 scope, also written up in README.md):
//   - No aromaticity. Every atom and bond is non-aromatic, so a perceived
//     component gets no rings: no plane-plane, atom-plane or group-plane
//     contacts, never typed aromatic. Ring C-C bonds (~1.39 A) read as
//     `Double`, so an unsubstituted aromatic CH carbon (two "double" ring
//     bonds = valence 4, not over-valent) gets no implicit H and isn't a
//     weak donor. Fixing that needs real Kekule/aromaticity perception.
//   - Neutral valences only, no charge model: an amine reads as neutral
//     (-NH2, three neighbours), so it's a donor and acceptor but not
//     pos-ionisable, unlike a CCD entry drawn as -NH3+ (e.g. lysine's NZ).
//     Real arpeggio protonates by pH via Open Babel.
//   - An atom bonded to another residue (a polymer link, a covalently
//     attached ligand) gets an implicit H in place of that bond, since
//     cross-residue bonds aren't perceived.
//   - Triple bonds read as `Double` (none of the fixtures have one).
//   - Valence repair only resolves a resonance group when its central atom
//     is over-valent; an amide whose carbonyl carbon has no third
//     substituent (formamide-like) keeps its C-N as `Double`.
//   - Open Babel's distance window is followed by a cleanup pass (delete
//     the longest bond on an atom over its element's *maximum* bond count,
//     or forming a <45 degree angle). Only the narrower bond-order repair
//     above is ported, so a badly distorted or clashing ligand can come
//     out over-bonded.
//   - One shared definition per comp_id, built from its first residue
//     instance. Atoms missing from that instance (disorder, partial
//     occupancy) are missing from every instance's typing.

use crate::config;
use crate::contacts::euclidean_distance;
use pdbtbx::{Atom, Residue, PDB};
use rspeggio_ccd::component::{BondOrder, CcdAtom, CcdBond, CcdComponent, ComponentType};
use std::collections::HashMap;

fn element_symbol(atom: &Atom) -> Option<String> {
    atom.element().map(|e| e.symbol().to_uppercase())
}

// Open Babel's own `OBMol::ConnectTheDots` bonding window (`mol.cpp`):
// bonded if 0.4 A < d < rcov_1 + rcov_2 + 0.45 A. The tolerance is needed,
// not cosmetic: real single bonds run slightly long of the bare radii sum
// (1UBQ's CA-CB/N-CA reach 1.028x it), so a bare-sum cutoff drops them.
const BOND_TOLERANCE: f64 = 0.45;
const MIN_BOND_DISTANCE: f64 = 0.4;

fn bond_order(element_1: &str, element_2: &str, distance: f64) -> Option<BondOrder> {
    let sum_cov_radii = config::covalent_radius(element_1)? + config::covalent_radius(element_2)?;
    if distance <= MIN_BOND_DISTANCE || distance >= sum_cov_radii + BOND_TOLERANCE {
        return None;
    }
    if distance < config::BOND_ORDER_DOUBLE_RATIO * sum_cov_radii {
        Some(BondOrder::Double)
    } else {
        Some(BondOrder::Single)
    }
}

// Standard neutral valence, for the elements repair and implicit-H
// counting cover. Sulfur is left out of repair (sulfonyl/sulfonate S is
// genuinely hypervalent) but still gets implicit H counted against 2, so
// a thiol gets its H and a sulfonyl S gets none.
fn standard_valence(element: &str) -> Option<u32> {
    match element {
        "C" => Some(4),
        "N" => Some(3),
        "O" | "S" => Some(2),
        _ => None,
    }
}

fn order_value(order: BondOrder) -> u32 {
    match order {
        BondOrder::Single | BondOrder::Aromatic => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
    }
}

struct PerceivedBond {
    i: usize,
    j: usize,
    order: BondOrder,
    distance: f64,
}

fn valence_sum(atom: usize, bonds: &[PerceivedBond]) -> u32 {
    bonds
        .iter()
        .filter(|b| b.i == atom || b.j == atom)
        .map(|b| order_value(b.order))
        .sum()
}

fn repair_over_valent(elements: &[String], bonds: &mut [PerceivedBond]) {
    loop {
        let demotion = (0..elements.len())
            .filter(|&a| elements[a] != "S")
            .filter_map(|a| Some((a, standard_valence(&elements[a])?)))
            .filter(|&(a, valence)| valence_sum(a, bonds) > valence)
            .find_map(|(a, _)| {
                bonds
                    .iter()
                    .enumerate()
                    .filter(|(_, b)| (b.i == a || b.j == a) && b.order == BondOrder::Double)
                    .max_by(|(_, x), (_, y)| x.distance.total_cmp(&y.distance))
                    .map(|(index, _)| index)
            });
        match demotion {
            Some(index) => bonds[index].order = BondOrder::Single,
            None => return,
        }
    }
}

// Builds a `CcdComponent` from `residue`'s real coordinates alone. Atoms
// are deduped by name (first wins), so alternate-location copies of the
// same atom don't become separate atoms bonded to each other. Atoms with
// no element can't be bonded and are skipped.
pub fn perceive_component(residue: &Residue) -> CcdComponent {
    let mut atoms: Vec<&Atom> = Vec::new();
    let mut elements: Vec<String> = Vec::new();
    for atom in residue.atoms() {
        if atoms.iter().any(|a| a.name() == atom.name()) {
            continue;
        }
        if let Some(element) = element_symbol(atom) {
            atoms.push(atom);
            elements.push(element);
        }
    }

    let mut bonds = Vec::new();
    for i in 0..atoms.len() {
        for j in (i + 1)..atoms.len() {
            let distance = euclidean_distance(atoms[i].pos(), atoms[j].pos());
            if let Some(order) = bond_order(&elements[i], &elements[j], distance) {
                bonds.push(PerceivedBond {
                    i,
                    j,
                    order,
                    distance,
                });
            }
        }
    }
    repair_over_valent(&elements, &mut bonds);

    let mut ccd_atoms: Vec<CcdAtom> = atoms
        .iter()
        .zip(&elements)
        .map(|(atom, element)| CcdAtom::new(atom.name().to_string(), element.clone(), false, false))
        .collect();
    let mut ccd_bonds: Vec<CcdBond> = bonds
        .iter()
        .map(|b| {
            CcdBond::new(
                atoms[b.i].name().to_string(),
                atoms[b.j].name().to_string(),
                b.order,
                false,
            )
        })
        .collect();

    for (a, atom) in atoms.iter().enumerate() {
        let Some(valence) = standard_valence(&elements[a]) else {
            continue;
        };
        for n in 1..=valence.saturating_sub(valence_sum(a, &bonds)) {
            let h_name = format!("{}_H{n}", atom.name());
            ccd_atoms.push(CcdAtom::new(h_name.clone(), "H".to_string(), false, false));
            ccd_bonds.push(CcdBond::new(
                atom.name().to_string(),
                h_name,
                BondOrder::Single,
                false,
            ));
        }
    }

    CcdComponent::new(ccd_atoms, ccd_bonds, ComponentType::BoundMolecule)
}

// Opt-in: for every comp_id in `pdb` with no entry in `components`,
// inserts a geometrically perceived one built from that comp_id's first
// residue instance. Entries already present (bundled or caller-loaded CCD
// files) are never replaced.
pub fn augment_with_geometric_fallback(pdb: &PDB, components: &mut HashMap<String, CcdComponent>) {
    for residue in pdb.residues() {
        let Some(comp_id) = residue.name() else {
            continue;
        };
        if !components.contains_key(comp_id) {
            components.insert(comp_id.to_string(), perceive_component(residue));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::typing::{bonded_neighbors, type_atom, AtomTypeBits, PHYSIOLOGICAL_PH};
    use rspeggio_ccd::common::common_components;
    use std::collections::BTreeSet;

    fn load(name: &str) -> PDB {
        let (pdb, _errors) = pdbtbx::open(format!("../../tests/fixtures/structures/{name}.cif"))
            .unwrap_or_else(|_| panic!("{name} should load"));
        pdb
    }

    fn first_residue<'a>(pdb: &'a PDB, comp_id: &str) -> &'a Residue {
        pdb.residues()
            .find(|r| r.name() == Some(comp_id))
            .unwrap_or_else(|| panic!("structure has at least one {comp_id}"))
    }

    fn atom<'a>(component: &'a CcdComponent, id: &str) -> &'a CcdAtom {
        component
            .atoms()
            .iter()
            .find(|a| a.atom_id() == id)
            .unwrap_or_else(|| panic!("{id} should be in the component"))
    }

    fn h_count(component: &CcdComponent, id: &str) -> usize {
        bonded_neighbors(atom(component, id), component)
            .iter()
            .filter(|(n, _)| n.element() == "H")
            .count()
    }

    fn heavy_bonds(component: &CcdComponent) -> Vec<(String, String, BondOrder)> {
        let is_h = |id: &str| atom(component, id).element() == "H";
        let mut bonds: Vec<_> = component
            .bonds()
            .iter()
            .filter(|b| !is_h(b.atom_id_1()) && !is_h(b.atom_id_2()))
            .map(|b| {
                let mut ids = [b.atom_id_1().to_string(), b.atom_id_2().to_string()];
                ids.sort();
                let [x, y] = ids;
                (x, y, *b.order())
            })
            .collect();
        bonds.sort_by(|a, b| (&a.0, &a.1).cmp(&(&b.0, &b.1)));
        bonds
    }

    #[test]
    fn a_real_alanines_perceived_bonds_and_hydrogens_match_its_real_ccd_entry() {
        // Ground truth: ALA's real (bundled) CCD entry, restricted to the
        // heavy atoms actually deposited in 1UBQ.
        let pdb = load("1UBQ");
        let ala = first_residue(&pdb, "ALA");
        let perceived = perceive_component(ala);

        let components = common_components();
        let real = &components["ALA"];
        let present: Vec<&str> = ala.atoms().map(|a| a.name()).collect();
        let real_heavy: Vec<_> = heavy_bonds(real)
            .into_iter()
            .filter(|(a, b, _)| present.contains(&a.as_str()) && present.contains(&b.as_str()))
            .collect();
        assert_eq!(heavy_bonds(&perceived), real_heavy);

        // Implicit H counts match the CCD's own (HA; HB1-3; H, H2).
        for id in ["CA", "CB", "N"] {
            assert_eq!(h_count(&perceived, id), h_count(real, id), "{id}");
        }
        // The carbonyl C really bonds to the next residue's N (or OXT in
        // the free monomer); with no cross-residue bonds perceived, that
        // slot gets an implicit H instead -- the documented polymer-link
        // limitation, asserted rather than hidden.
        assert_eq!(h_count(&perceived, "C"), 1);
        assert_eq!(h_count(real, "C"), 0);

        assert_eq!(perceived.component_type(), ComponentType::BoundMolecule);
        assert!(perceived.atoms().iter().all(|a| !a.aromatic()));
    }

    #[test]
    fn a_real_amide_c_n_bond_reads_as_double_before_valence_repair() {
        // MET1 C -> GLN2 N in 1UBQ: a real peptide bond (single in the
        // CCD, resonance-shortened) at 1.299 A, ratio 0.884 to the summed
        // C+N radii -- inside the C=O double-bond range, so the initial
        // distance pass alone calls it double. Valence repair is what
        // resolves this (see the ASN test below).
        let pdb = load("1UBQ");
        let met1 = pdb.residues().find(|r| r.id().0 == 1).expect("MET1");
        let gln2 = pdb.residues().find(|r| r.id().0 == 2).expect("GLN2");
        let c = met1.atoms().find(|a| a.name() == "C").expect("MET1 C");
        let n = gln2.atoms().find(|a| a.name() == "N").expect("GLN2 N");

        let distance = euclidean_distance(c.pos(), n.pos());
        assert!((distance - 1.299).abs() < 0.001);
        assert_eq!(bond_order("C", "N", distance), Some(BondOrder::Double));
    }

    #[test]
    fn valence_repair_recovers_a_real_side_chain_amide() {
        // ASN's CG-ND2 reads as double from distance alone (CG would then
        // carry 5 bond orders); repair demotes the longer C-N, leaving
        // CG=OD1 as the double bond, so amide perception finds the same
        // group it finds in ASN's real CCD entry.
        let pdb = load("1UBQ");
        let perceived = perceive_component(first_residue(&pdb, "ASN"));

        let groups = crate::rings::perceive_amide_groups(&perceived);
        assert_eq!(
            groups,
            vec![crate::rings::AmideGroup {
                nitrogen_id: "ND2".to_string(),
                carbon_id: "CG".to_string(),
                oxygen_id: "OD1".to_string(),
                other_carbon_id: "CB".to_string(),
            }]
        );
        assert_eq!(h_count(&perceived, "ND2"), 2);
    }

    #[test]
    fn valence_repair_recovers_a_real_carboxylate() {
        // Both of a real deprotonated carboxylate's C-O bonds read as
        // double (114 of 115 across every fixture here -- resonance makes
        // them equal, ~1.25 A). Repair keeps the shorter one double, the
        // single-bonded O gets an implicit H (the CCD's own free-acid
        // drawing), and typing's pH rule deprotonates it as usual.
        let pdb = load("1UBQ");
        let asp = first_residue(&pdb, "ASP");
        let perceived = perceive_component(asp);

        let cg_bonds: Vec<_> = bonded_neighbors(atom(&perceived, "CG"), &perceived)
            .into_iter()
            .filter(|(n, _)| n.element() == "O")
            .map(|(n, order)| (n.atom_id().to_string(), *order))
            .collect();
        let doubles: Vec<_> = cg_bonds
            .iter()
            .filter(|(_, o)| *o == BondOrder::Double)
            .collect();
        assert_eq!(doubles.len(), 1, "{cg_bonds:?}");

        let pos = |id: &str| asp.atoms().find(|a| a.name() == id).unwrap().pos();
        let cg = pos("CG");
        let shorter = ["OD1", "OD2"]
            .into_iter()
            .min_by(|a, b| {
                euclidean_distance(cg, pos(a)).total_cmp(&euclidean_distance(cg, pos(b)))
            })
            .unwrap();
        assert_eq!(doubles[0].0, shorter);

        for id in ["OD1", "OD2"] {
            let bits = type_atom(atom(&perceived, id), &perceived, PHYSIOLOGICAL_PH);
            assert!(bits.contains(AtomTypeBits::NEG_IONISABLE), "{id}");
            assert!(bits.contains(AtomTypeBits::HBOND_ACCEPTOR), "{id}");
            assert!(!bits.contains(AtomTypeBits::HBOND_DONOR), "{id} at pH 7.4");
        }
    }

    #[test]
    fn non_bonded_atoms_in_the_same_residue_get_no_bond() {
        // ALA's N and C are two bonds apart (via CA), ~2.4-2.5 A -- well
        // beyond the 1.92 A cutoff (1.47 A summed radii + 0.45 A).
        let pdb = load("1UBQ");
        let perceived = perceive_component(first_residue(&pdb, "ALA"));
        let bonded = |a: &str, b: &str| {
            bonded_neighbors(atom(&perceived, a), &perceived)
                .iter()
                .any(|(n, _)| n.atom_id() == b)
        };
        assert!(!bonded("N", "C"));
        assert!(bonded("N", "CA"));
    }

    #[test]
    fn perceived_typing_matches_real_ccd_typing_except_for_known_limitations() {
        // Every amino-acid type in 1CA2 (first instance of each), treated
        // as if unknown: perceived typing and H counts vs. the real CCD
        // entry, for side-chain atoms plus the carbonyl O (backbone N/C
        // are polymer-link atoms, see the ALA test). 77 of 127 atoms
        // match exactly; every mismatch is pinned here by cause.
        let pdb = load("1CA2");
        let components = common_components();

        let mut expected: BTreeSet<String> = BTreeSet::new();
        // The CCD's free-monomer OXT makes every residue's own O look like
        // a carboxylate oxygen (NEG_IONISABLE); perceived, with no OXT in
        // an internal residue, doesn't. A CCD-side quirk, not this module's.
        for res in [
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS",
            "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        ] {
            expected.insert(format!("{res} O"));
        }
        // No aromaticity in v1: ring atoms lose AROMATIC, ring N-H lose
        // donor status, some ring CH lose their implicit H.
        for id in [
            "HIS CG", "HIS ND1", "HIS CD2", "HIS CE1", "HIS NE2", "TRP CG", "TRP CD1", "TRP CD2",
            "TRP NE1", "TRP CE2", "TRP CE3", "TRP CZ2", "TRP CZ3", "TRP CH2", "TYR CG", "TYR CD1",
            "TYR CD2", "TYR CE1", "TYR CE2", "TYR CZ", "PHE CG", "PHE CD1", "PHE CD2", "PHE CE1",
            "PHE CE2", "PHE CZ",
        ] {
            expected.insert(id.to_string());
        }
        // No charge model: neutral -NH2 (donor + acceptor), not -NH3+.
        expected.insert("LYS NZ".to_string());
        // Resonance-equivalent carboxylate O's: the other one is drawn as
        // C=O (both still NEG_IONISABLE acceptors).
        expected.insert("GLU OE1".to_string());
        expected.insert("GLU OE2".to_string());
        // Guanidinium: typing identical, but NH2 keeps the double bond so
        // gets 1 implicit H instead of the CCD's charged form's 2.
        expected.insert("ARG NH2".to_string());

        let mut seen = BTreeSet::new();
        let mut mismatched = BTreeSet::new();
        let mut total = 0;
        for residue in pdb.residues() {
            let name = residue.name().unwrap_or_default();
            if !crate::config::STANDARD_AMINO_ACIDS.contains(&name) || !seen.insert(name) {
                continue;
            }
            let real = &components[name];
            let perceived = perceive_component(residue);
            for structure_atom in residue.atoms() {
                let id = structure_atom.name();
                if ["N", "C", "OXT"].contains(&id) {
                    continue;
                }
                let Some(real_atom) = real.atoms().iter().find(|a| a.atom_id() == id) else {
                    continue;
                };
                total += 1;
                let same_type = type_atom(real_atom, real, PHYSIOLOGICAL_PH)
                    == type_atom(atom(&perceived, id), &perceived, PHYSIOLOGICAL_PH);
                if !same_type || h_count(real, id) != h_count(&perceived, id) {
                    mismatched.insert(format!("{name} {id}"));
                }
            }
        }

        assert_eq!(seen.len(), 20, "1CA2 has all 20 standard amino acids");
        assert_eq!(total, 127);
        assert_eq!(mismatched, expected);
    }

    #[test]
    fn augmenting_fills_in_a_missing_component_and_keeps_existing_ones() {
        let pdb = load("1UBQ");
        let mut components = common_components();
        components.remove("ALA");
        let gly_bonds_before = components["GLY"].bonds().len();

        augment_with_geometric_fallback(&pdb, &mut components);

        let ala = components.get("ALA").expect("ALA should now be perceived");
        assert_eq!(ala.component_type(), ComponentType::BoundMolecule);
        assert_eq!(
            components["GLY"].bonds().len(),
            gly_bonds_before,
            "a real CCD entry must never be replaced by a perceived one"
        );
        assert_eq!(
            components["GLY"].component_type(),
            ComponentType::Polypeptide
        );
    }

    #[test]
    fn without_augmenting_an_unknown_component_still_fails_loudly() {
        let pdb = load("1UBQ");
        let mut components = common_components();
        components.remove("ALA");
        let selection = crate::selection::SelectionContext::whole_structure(&pdb);

        let result = crate::export::export_atom_atom_contacts(
            &pdb,
            &components,
            &selection,
            PHYSIOLOGICAL_PH,
        );
        assert!(result.is_err());
    }

    fn ala_contacts(perceived: bool) -> Vec<(String, Vec<&'static str>, &'static str)> {
        let pdb = load("1UBQ");
        let mut components = common_components();
        if perceived {
            components.remove("ALA");
            augment_with_geometric_fallback(&pdb, &mut components);
        }
        let selection = crate::selection::SelectionContext::whole_structure(&pdb);
        let entries = crate::export::export_atom_atom_contacts(
            &pdb,
            &components,
            &selection,
            PHYSIOLOGICAL_PH,
        )
        .expect("every component is known (bundled or perceived)");

        let mut contacts: Vec<_> = entries
            .into_iter()
            .filter_map(|e| {
                let ala = [&e.bgn, &e.end]
                    .into_iter()
                    .find(|a| a.label_comp_id == "ALA")?
                    .label_comp_type;
                let key = format!(
                    "{}{}-{}{}",
                    e.bgn.auth_seq_id, e.bgn.auth_atom_id, e.end.auth_seq_id, e.end.auth_atom_id
                );
                Some((key, e.contact, ala))
            })
            .collect();
        contacts.sort();
        contacts
    }

    #[test]
    fn an_augmented_unknown_component_exports_nearly_the_same_contacts_as_its_real_ccd_entry() {
        // 1UBQ with ALA's CCD entry removed and perceived instead, vs. the
        // real entry: the same 76 ALA contacts, 73 labelled identically
        // (all 9 hbond, 11 polar, 3 weak_hbond). The 3 that differ are
        // both documented limitations showing up on real data.
        let real = ala_contacts(false);
        let perceived = ala_contacts(true);

        assert_eq!(real.len(), 76);
        let keys =
            |c: &[(String, Vec<&str>, &str)]| c.iter().map(|x| x.0.clone()).collect::<Vec<_>>();
        assert_eq!(keys(&perceived), keys(&real));
        assert!(real.iter().all(|c| c.2 == "P"));
        assert!(perceived.iter().all(|c| c.2 == "B"));

        let differing: Vec<_> = real
            .iter()
            .zip(&perceived)
            .filter(|(r, p)| r.1 != p.1)
            .map(|(r, p)| (r.0.as_str(), r.1.clone(), p.1.clone()))
            .collect();
        assert_eq!(
            differing,
            vec![
                // The phantom implicit H on ALA's backbone C (it really
                // bonds to the next residue's N) makes it a weak donor.
                (
                    "28C-30N",
                    vec!["vdw_clash"],
                    vec!["vdw_clash", "weak_polar"]
                ),
                (
                    "46C-48N",
                    vec!["vdw_clash"],
                    vec!["vdw_clash", "weak_polar"]
                ),
                // The real CCD's free-monomer OXT makes ALA46's internal O
                // look like a carboxylate; perceived correctly doesn't.
                (
                    "46O-48NZ",
                    vec!["proximal", "hbond", "ionic", "polar"],
                    vec!["proximal", "hbond", "polar"]
                ),
            ]
        );
    }
}
