// crates/rspeggio-core/src/hydrogenate.rs
//
// Analytically places a hydrogen the CCD says exists but a real structure
// doesn't resolve (the normal case -- most X-ray structures have no
// hydrogens at all). This is the same problem OpenBabel's
// `AddHydrogens`/`OBBuilder::GetNewBondVector` solves, but by the same
// local, analytic vector geometry OpenBabel itself uses (bond length from
// covalent radii, direction from the atom's existing real neighbor
// positions + standard bond angles for its hybridization) -- not by
// linking OpenBabel. See the `hydrogenation_diagram` artifact for the
// overall decision flow; this module builds it case by case.
//
// A systematic sweep of every bundled/fixture CCD file's real bond graphs
// (comp_id, neighbor_count, hydrogen_count) found: 1-neighbor/1H (63
// occurrences), 1-neighbor/2H planar (30), 1-neighbor/3H methyl (11),
// 2-neighbor/2H sp3 CH2 (42), 2-neighbor/1H sp2 (45), 3-neighbor/1H (56)
// -- all covered below -- versus 0-neighbor/2H water (1 occurrence, and
// deliberately *not* handled here: see below).
//
// Water is the wrong problem for this module, not a gap in it: its
// oxygen has no heavy-atom neighbor at all to anchor a bond direction
// from, so placing its two hydrogens needs the orientation of the whole
// molecule -- which depends on the surrounding hydrogen-bond network, not
// the bond graph. Real pdbe-arpeggio doesn't solve this either --
// `interactions.py:791` explicitly special-cases water to skip the
// H-based check entirely, treating its oxygen as unconditionally
// donor-and-acceptor-capable. The fix belongs in the future
// feature-contact code as that same kind of special case, not here.
//
// Not yet implemented: the sp3 2-known-neighbor case with only *one*
// missing hydrogen (a secondary CH, distinct from the CH2 case below,
// which needs two); a pyramidal (non-planar) 1-neighbor 2-hydrogen case,
// for a genuinely sp3 primary amine rather than a resonance-planar amide
// (real side-chain amides like Asn/Gln, and adenine's exocyclic N6, are
// both planar).
//
// Every function below does the geometry once given real 3D positions --
// none of them walk the structure themselves to find a reference atom
// (e.g. the previous residue's real C, via sequence adjacency, for the
// trans-to-reference case). That lookup is still a caller's job.

type Point = (f64, f64, f64);

pub enum Hybridization {
    Sp2,
    Sp3,
}

fn subtract(a: Point, b: Point) -> Point {
    (a.0 - b.0, a.1 - b.1, a.2 - b.2)
}

fn length(v: Point) -> f64 {
    (v.0 * v.0 + v.1 * v.1 + v.2 * v.2).sqrt()
}

// The angle a-vertex-b, in degrees. Used both by this module's own tests
// and, via `hydrogenate::angle_degrees`, by `features.rs`'s real hbond
// angle check.
pub fn angle_degrees(a: Point, vertex: Point, b: Point) -> f64 {
    let v1 = normalize(subtract(a, vertex));
    let v2 = normalize(subtract(b, vertex));
    dot(v1, v2).clamp(-1.0, 1.0).acos().to_degrees()
}

fn normalize(v: Point) -> Point {
    let len = length(v);
    (v.0 / len, v.1 / len, v.2 / len)
}

fn dot(a: Point, b: Point) -> f64 {
    a.0 * b.0 + a.1 * b.1 + a.2 * b.2
}

fn cross(a: Point, b: Point) -> Point {
    (
        a.1 * b.2 - a.2 * b.1,
        a.2 * b.0 - a.0 * b.2,
        a.0 * b.1 - a.1 * b.0,
    )
}

// The component of `v` perpendicular to `axis_unit` (which must already be
// a unit vector), normalized. Ill-conditioned when `v` is (nearly)
// parallel to `axis_unit` -- callers are responsible for picking a `v`
// that isn't.
fn perpendicular_component(v: Point, axis_unit: Point) -> Point {
    let along = dot(v, axis_unit);
    normalize((
        v.0 - axis_unit.0 * along,
        v.1 - axis_unit.1 * along,
        v.2 - axis_unit.2 * along,
    ))
}

// Places the fourth substituent on a tetrahedral (sp3) center, given the
// real 3D positions of the center and its other three neighbors. This
// case is fully determined by geometry alone: an ideal sp3 center's four
// bond vectors sum to (approximately) zero, so the one missing direction
// is just the negative sum of the three known unit bond vectors -- no
// hybridization choice or reference atom needed, unlike the 1- and
// 2-neighbor cases.
pub fn place_fourth_tetrahedral_substituent(
    center: Point,
    neighbor_1: Point,
    neighbor_2: Point,
    neighbor_3: Point,
    bond_length: f64,
) -> Point {
    let u1 = normalize(subtract(neighbor_1, center));
    let u2 = normalize(subtract(neighbor_2, center));
    let u3 = normalize(subtract(neighbor_3, center));

    let sum = (u1.0 + u2.0 + u3.0, u1.1 + u2.1 + u3.1, u1.2 + u2.2 + u3.2);
    let direction = normalize((-sum.0, -sum.1, -sum.2));

    (
        center.0 + direction.0 * bond_length,
        center.1 + direction.1 * bond_length,
        center.2 + direction.2 * bond_length,
    )
}

// Places a substituent given two known heavy neighbors, for a planar
// (sp2) center -- e.g. an aromatic ring nitrogen's hydrogen, with both
// ring-carbon neighbors already resolved (HIS's ND1, bonded to CG and
// CE1). Fully determined, same construction as the three-neighbor case:
// the missing direction is the in-plane bisector of the two known bonds,
// negated -- pointing away from both, in the plane they define. Only
// valid for sp2 centers; an sp3 center with two known neighbors instead
// needs the direction perpendicular to that plane (not yet implemented).
pub fn place_planar_two_neighbor_substituent(
    center: Point,
    neighbor_1: Point,
    neighbor_2: Point,
    bond_length: f64,
) -> Point {
    let u1 = normalize(subtract(neighbor_1, center));
    let u2 = normalize(subtract(neighbor_2, center));

    let sum = (u1.0 + u2.0, u1.1 + u2.1, u1.2 + u2.2);
    let direction = normalize((-sum.0, -sum.1, -sum.2));

    (
        center.0 + direction.0 * bond_length,
        center.1 + direction.1 * bond_length,
        center.2 + direction.2 * bond_length,
    )
}

// Places a substituent given exactly one known heavy neighbor, at the
// ideal bond angle for the given hybridization (109.47 degrees for sp3,
// 120 for sp2). Underdetermined by geometry alone: with only one real
// constraint, there's a whole cone of directions at the correct angle,
// and this picks an arbitrary but deterministic one on that cone (via an
// arbitrary reference axis), rather than using a second-shell atom to fix
// the azimuth the way OpenBabel does for e.g. keeping an amide N-H trans
// to its carbonyl. That's an acceptable simplification for genuinely
// freely-rotating single substituents (a hydroxyl's H has real
// rotational freedom too); it's the wrong tool for cases where the
// specific azimuth is chemically fixed, which isn't implemented yet.
pub fn place_arbitrary_single_neighbor_substituent(
    center: Point,
    known_neighbor: Point,
    hybridization: Hybridization,
    bond_length: f64,
) -> Point {
    let bond1 = normalize(subtract(known_neighbor, center));
    let target_angle = match hybridization {
        Hybridization::Sp2 => 120.0_f64.to_radians(),
        Hybridization::Sp3 => 109.471_22_f64.to_radians(),
    };

    // An arbitrary reference axis, picked so it's never (nearly) parallel
    // to bond1 -- otherwise the perpendicular component below would be
    // ill-conditioned (close to the zero vector).
    let axis = if dot(bond1, (0.0, 0.0, 1.0)).abs() < 0.9 {
        (0.0, 0.0, 1.0)
    } else {
        (1.0, 0.0, 0.0)
    };
    let perpendicular = perpendicular_component(axis, bond1);

    let (sin_a, cos_a) = target_angle.sin_cos();
    let direction = (
        bond1.0 * cos_a + perpendicular.0 * sin_a,
        bond1.1 * cos_a + perpendicular.1 * sin_a,
        bond1.2 * cos_a + perpendicular.2 * sin_a,
    );

    (
        center.0 + direction.0 * bond_length,
        center.1 + direction.1 * bond_length,
        center.2 + direction.2 * bond_length,
    )
}

// Places a substituent given exactly one known heavy neighbor plus a
// second-shell reference atom (bonded to the known neighbor, not to
// `center` itself) -- e.g. a backbone amide N-H: `center` is N,
// `known_neighbor` is CA, and `reference` is the *previous* residue's
// real carbonyl carbon, the actual peptide-bond partner that defines the
// amide plane. Same angle construction as the arbitrary-azimuth version,
// but the azimuth is fixed by placing the new substituent trans (opposite
// side, viewed along the known_neighbor-center bond) to the reference,
// rather than picked arbitrarily -- this is the OpenBabel
// `GetNewBondVector` behavior its own comments describe as "make sure to
// place the new atom trans to a-2".
pub fn place_single_neighbor_substituent_trans_to_reference(
    center: Point,
    known_neighbor: Point,
    reference: Point,
    hybridization: Hybridization,
    bond_length: f64,
) -> Point {
    let bond1 = normalize(subtract(known_neighbor, center));
    let target_angle = match hybridization {
        Hybridization::Sp2 => 120.0_f64.to_radians(),
        Hybridization::Sp3 => 109.471_22_f64.to_radians(),
    };

    // The reference's azimuthal direction, as seen from known_neighbor
    // (the pivot the bond1 axis passes through) -- then place the new
    // substituent's azimuth on the opposite side (negate it).
    let toward_reference = subtract(reference, known_neighbor);
    let perpendicular_toward_reference = perpendicular_component(toward_reference, bond1);
    let perpendicular = (
        -perpendicular_toward_reference.0,
        -perpendicular_toward_reference.1,
        -perpendicular_toward_reference.2,
    );

    let (sin_a, cos_a) = target_angle.sin_cos();
    let direction = (
        bond1.0 * cos_a + perpendicular.0 * sin_a,
        bond1.1 * cos_a + perpendicular.1 * sin_a,
        bond1.2 * cos_a + perpendicular.2 * sin_a,
    );

    (
        center.0 + direction.0 * bond_length,
        center.1 + direction.1 * bond_length,
        center.2 + direction.2 * bond_length,
    )
}

// Places both hydrogens of a planar (resonance-delocalized) primary
// amide amine, e.g. asparagine's ND2 or glutamine's NE2: `center` is the
// nitrogen, `known_neighbor` is its one heavy neighbor (CG/CD), and
// `reference` is that neighbor's other real substituent (e.g. the
// amide's OD1/OE1, part of the same conjugated C=O/C-N plane) used to
// orient the plane -- analogous to the trans-to-reference case, but here
// the reference fixes the shared plane both hydrogens lie in, rather than
// picking a single opposite side.
//
// All three substituents on a trigonal planar center (known_neighbor, H1,
// H2) sit in one plane, 120 degrees apart from each other in it. With
// `bond1` (toward known_neighbor) at in-plane angle 0, the other two are
// placed at +120 and -120 degrees within that same plane -- which by
// construction also puts them at exactly 120 degrees from bond1 and from
// each other.
pub fn place_planar_amine_hydrogens(
    center: Point,
    known_neighbor: Point,
    reference: Point,
    bond_length: f64,
) -> (Point, Point) {
    let bond1 = normalize(subtract(known_neighbor, center));
    let toward_reference = subtract(reference, known_neighbor);
    let perp = perpendicular_component(toward_reference, bond1);

    let (sin_a, cos_a) = 120.0_f64.to_radians().sin_cos();
    let place = |perp_sign: f64| {
        let direction = (
            bond1.0 * cos_a + perp.0 * sin_a * perp_sign,
            bond1.1 * cos_a + perp.1 * sin_a * perp_sign,
            bond1.2 * cos_a + perp.2 * sin_a * perp_sign,
        );
        (
            center.0 + direction.0 * bond_length,
            center.1 + direction.1 * bond_length,
            center.2 + direction.2 * bond_length,
        )
    };

    (place(1.0), place(-1.0))
}

// Places both hydrogens of a tetrahedral (sp3) center with exactly two
// known heavy neighbors -- a methylene (CH2) group, e.g. serine's CB
// (bonded to CA and OG) or glycerol's terminal carbons. The confirmed
// highest-volume real gap this module had: 42 real occurrences across
// the bundled common set and fixture ligands, more than any other single
// pattern -- CH2 is one of the most ordinary motifs in organic chemistry,
// not a corner case.
//
// Standard tetrahedral construction: with `bisector` the (normalized) sum
// of the two known bond directions and `perp` perpendicular to the plane
// they define, the two new substituents are placed symmetrically off
// *negative* bisector (away from both known neighbors), split into the
// +perp/-perp directions by half the tetrahedral angle. This guarantees
// the angle between the two placed hydrogens is exactly the tetrahedral
// angle regardless of how close the two *real* input neighbors are to
// ideal geometry themselves (see the real-data test) -- the same kind of
// invariant `place_fourth_tetrahedral_substituent` guarantees for its own
// construction.
pub fn place_tetrahedral_pair_from_two_neighbors(
    center: Point,
    neighbor_1: Point,
    neighbor_2: Point,
    bond_length: f64,
) -> (Point, Point) {
    let u1 = normalize(subtract(neighbor_1, center));
    let u2 = normalize(subtract(neighbor_2, center));

    let bisector = normalize((u1.0 + u2.0, u1.1 + u2.1, u1.2 + u2.2));
    let perp = normalize(cross(u1, u2));

    let half_tetrahedral_angle = (109.471_22_f64 / 2.0).to_radians();
    let (sin_a, cos_a) = half_tetrahedral_angle.sin_cos();

    let place = |perp_sign: f64| {
        let direction = (
            -bisector.0 * cos_a + perp.0 * sin_a * perp_sign,
            -bisector.1 * cos_a + perp.1 * sin_a * perp_sign,
            -bisector.2 * cos_a + perp.2 * sin_a * perp_sign,
        );
        (
            center.0 + direction.0 * bond_length,
            center.1 + direction.1 * bond_length,
            center.2 + direction.2 * bond_length,
        )
    };

    (place(1.0), place(-1.0))
}

// Places all three hydrogens of a methyl (CH3) group: a tetrahedral (sp3)
// center with exactly one known heavy neighbor. Underdetermined in the
// same rotational sense as `place_arbitrary_single_neighbor_substituent`
// -- a real methyl genuinely rotates freely, so there's no chemically
// "correct" azimuth to recover, only a consistent one. Picks an arbitrary
// reference axis the same way that function does, then places the three
// hydrogens 120 degrees apart in azimuth around the bond axis, each at
// the tetrahedral polar angle from it.
pub fn place_methyl_hydrogens(
    center: Point,
    known_neighbor: Point,
    bond_length: f64,
) -> (Point, Point, Point) {
    let bond1 = normalize(subtract(known_neighbor, center));

    let axis = if dot(bond1, (0.0, 0.0, 1.0)).abs() < 0.9 {
        (0.0, 0.0, 1.0)
    } else {
        (1.0, 0.0, 0.0)
    };
    let perp1 = perpendicular_component(axis, bond1);
    let perp2 = cross(bond1, perp1); // already unit length: bond1 and perp1 are orthonormal

    let polar = 109.471_22_f64.to_radians();
    let (sin_p, cos_p) = polar.sin_cos();

    let place = |azimuth_deg: f64| {
        let (sin_az, cos_az) = azimuth_deg.to_radians().sin_cos();
        let direction = (
            bond1.0 * cos_p + (perp1.0 * cos_az + perp2.0 * sin_az) * sin_p,
            bond1.1 * cos_p + (perp1.1 * cos_az + perp2.1 * sin_az) * sin_p,
            bond1.2 * cos_p + (perp1.2 * cos_az + perp2.2 * sin_az) * sin_p,
        );
        (
            center.0 + direction.0 * bond_length,
            center.1 + direction.1 * bond_length,
            center.2 + direction.2 * bond_length,
        )
    };

    (place(0.0), place(120.0), place(240.0))
}

#[cfg(test)]
mod tests {
    use super::*;
    use pdbtbx::{ContainsAtomConformer, ContainsAtomConformerResidue};

    #[test]
    fn places_alanines_ha_from_real_1ubq_geometry_at_the_right_length_and_angles() {
        // Real MET1 heavy-atom positions, read directly from
        // tests/fixtures/structures/1UBQ.cif -- an ordinary X-ray structure
        // with no hydrogens at all, which is exactly the situation this
        // function exists for.
        let (pdb, _errors) = pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif")
            .expect("1UBQ fixture should parse");
        let met1 = |name: &str| -> Point {
            pdb.atoms_with_hierarchy()
                .find(|h| h.residue().id().0 == 1 && h.atom().name() == name)
                .unwrap_or_else(|| panic!("MET1 {name} should be present"))
                .atom()
                .pos()
        };
        let n = met1("N");
        let ca = met1("CA");
        let c = met1("C");
        let cb = met1("CB");

        // Real MET1 backbone angles at CA are themselves non-ideal --
        // N-CA-C is ~107.0 degrees, N-CA-CB ~113.9, C-CA-CB ~111.2 (all
        // confirmed directly from these same coordinates), a natural ~7
        // degree spread around the textbook 109.47. So the placed
        // direction won't sit at a uniform tetrahedral angle from all
        // three real neighbors either -- that would only hold for a
        // perfectly symmetric input. What the construction *does*
        // guarantee, by definition, is that the placed direction is
        // exactly antiparallel to the sum of the three real unit bond
        // vectors -- not that all four sum to zero (that only holds when
        // the sum happens to have unit magnitude, true only for
        // perfectly symmetric input).
        let bond_length = 1.09; // typical C-H bond length
        let ha = place_fourth_tetrahedral_substituent(ca, n, c, cb, bond_length);

        let actual_length = length(subtract(ha, ca));
        assert!(
            (actual_length - bond_length).abs() < 1e-9,
            "bond length should be exact by construction, got {actual_length}"
        );

        let u_ha = normalize(subtract(ha, ca));
        let u_n = normalize(subtract(n, ca));
        let u_c = normalize(subtract(c, ca));
        let u_cb = normalize(subtract(cb, ca));
        let sum = (
            u_n.0 + u_c.0 + u_cb.0,
            u_n.1 + u_c.1 + u_cb.1,
            u_n.2 + u_c.2 + u_cb.2,
        );
        let cos_angle_to_sum = dot(u_ha, normalize(sum));
        assert!(
            (cos_angle_to_sum - (-1.0)).abs() < 1e-9,
            "placed direction should be exactly antiparallel to the sum of the three real \
             unit bond vectors, got cos(angle) = {cos_angle_to_sum}"
        );

        // Still a real, independent sanity bound: the placed direction
        // shouldn't be wildly wrong (e.g. pointing back through an
        // existing neighbor).
        for neighbor in [n, c, cb] {
            let angle = angle_degrees(ha, ca, neighbor);
            assert!(
                (80.0..140.0).contains(&angle),
                "HA-CA-neighbor angle {angle} is nowhere near a plausible tetrahedral range"
            );
        }
    }

    #[test]
    fn a_symmetric_center_places_the_fourth_substituent_straight_up() {
        // Three neighbors arranged symmetrically in the xy-plane below the
        // center -- the missing 4th direction should point straight up
        // (+Z), a case simple enough to check by hand rather than just by
        // its own invariants.
        let center = (0.0, 0.0, 0.0);
        let r = 1.0;
        let half_sqrt_3 = 3.0_f64.sqrt() / 2.0;
        let n1 = (r, 0.0, -0.5);
        let n2 = (-r * 0.5, r * half_sqrt_3, -0.5);
        let n3 = (-r * 0.5, -r * half_sqrt_3, -0.5);

        let placed = place_fourth_tetrahedral_substituent(center, n1, n2, n3, 1.0);

        assert!(placed.0.abs() < 1e-9, "x should cancel out: {placed:?}");
        assert!(placed.1.abs() < 1e-9, "y should cancel out: {placed:?}");
        assert!(
            placed.2 > 0.0,
            "should point away from the three neighbors: {placed:?}"
        );
    }

    #[test]
    fn places_histidines_nd1_hydrogen_from_real_1ubq_ring_geometry() {
        // Real HIS68 heavy-atom positions from 1UBQ. ND1 is bonded to CG
        // (single) and CE1 (double, ring), both already resolved -- HD1
        // is the missing hydrogen this places.
        let (pdb, _errors) = pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif")
            .expect("1UBQ fixture should parse");
        let his68 = |name: &str| -> Point {
            pdb.atoms_with_hierarchy()
                .find(|h| {
                    h.residue().id().0 == 68
                        && h.residue().name() == Some("HIS")
                        && h.atom().name() == name
                })
                .unwrap_or_else(|| panic!("HIS68 {name} should be present"))
                .atom()
                .pos()
        };
        let nd1 = his68("ND1");
        let cg = his68("CG");
        let ce1 = his68("CE1");

        let bond_length = 1.02; // N-H: covalent radii 0.71 + 0.31
        let hd1 = place_planar_two_neighbor_substituent(nd1, cg, ce1, bond_length);

        let actual_length = length(subtract(hd1, nd1));
        assert!(
            (actual_length - bond_length).abs() < 1e-9,
            "bond length should be exact by construction, got {actual_length}"
        );

        // Fully determined by construction (unlike the 3-neighbor case,
        // this uses only 2 real constraints, so both angles land exactly
        // where the in-plane-bisector formula puts them) -- verify the
        // direction is genuinely opposite both known bonds, not just
        // "somewhere plausible".
        for neighbor in [cg, ce1] {
            let angle = angle_degrees(hd1, nd1, neighbor);
            assert!(
                angle > 90.0,
                "HD1 should point away from {neighbor:?}, got angle {angle}"
            );
        }

        // The placed hydrogen should lie in the same plane as ND1's two
        // real neighbors (a genuine planar/aromatic ring), i.e. the
        // normal to (CG, ND1, CE1) should be ~perpendicular to ND1-HD1.
        let normal = {
            let a = subtract(cg, nd1);
            let b = subtract(ce1, nd1);
            (
                a.1 * b.2 - a.2 * b.1,
                a.2 * b.0 - a.0 * b.2,
                a.0 * b.1 - a.1 * b.0,
            )
        };
        let cos_to_normal = dot(normalize(subtract(hd1, nd1)), normalize(normal));
        assert!(
            cos_to_normal.abs() < 1e-6,
            "HD1 should lie in the ring plane, got cos(angle to normal) = {cos_to_normal}"
        );
    }

    #[test]
    fn places_serines_og_hydrogen_from_real_1ubq_geometry() {
        // Real SER20 heavy-atom positions from 1UBQ. OG is bonded only to
        // CB -- HG is the missing hydroxyl hydrogen this places.
        let (pdb, _errors) = pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif")
            .expect("1UBQ fixture should parse");
        let ser20 = |name: &str| -> Point {
            pdb.atoms_with_hierarchy()
                .find(|h| {
                    h.residue().id().0 == 20
                        && h.residue().name() == Some("SER")
                        && h.atom().name() == name
                })
                .unwrap_or_else(|| panic!("SER20 {name} should be present"))
                .atom()
                .pos()
        };
        let og = ser20("OG");
        let cb = ser20("CB");

        let bond_length = 0.97; // O-H: covalent radii 0.66 + 0.31
        let hg =
            place_arbitrary_single_neighbor_substituent(og, cb, Hybridization::Sp3, bond_length);

        let actual_length = length(subtract(hg, og));
        assert!(
            (actual_length - bond_length).abs() < 1e-9,
            "bond length should be exact by construction, got {actual_length}"
        );

        // Fully determined here too (the formula fixes the angle exactly
        // at construction time, even though the azimuth around CB-OG is
        // arbitrary).
        let angle = angle_degrees(hg, og, cb);
        assert!(
            (angle - 109.471_22).abs() < 1e-6,
            "OG-HG should sit at the exact sp3 angle from CB by construction, got {angle}"
        );
    }

    #[test]
    fn places_a_backbone_amide_nh_trans_to_the_real_previous_residues_carbonyl() {
        // Real cross-residue peptide-bond geometry from 1UBQ: GLN2's N
        // (bonded only to CA within its own residue) and MET1's real C --
        // the actual peptide-bond partner that defines the amide plane,
        // not an intra-residue stand-in.
        let (pdb, _errors) = pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif")
            .expect("1UBQ fixture should parse");
        let atom_pos = |seq_id: isize, comp_id: &str, name: &str| -> Point {
            pdb.atoms_with_hierarchy()
                .find(|h| {
                    h.residue().id().0 == seq_id
                        && h.residue().name() == Some(comp_id)
                        && h.atom().name() == name
                })
                .unwrap_or_else(|| panic!("{comp_id}{seq_id} {name} should be present"))
                .atom()
                .pos()
        };
        let n = atom_pos(2, "GLN", "N");
        let ca = atom_pos(2, "GLN", "CA");
        let prev_c = atom_pos(1, "MET", "C");

        let bond_length = 1.02; // N-H: covalent radii 0.71 + 0.31
        let h = place_single_neighbor_substituent_trans_to_reference(
            n,
            ca,
            prev_c,
            Hybridization::Sp2,
            bond_length,
        );

        let actual_length = length(subtract(h, n));
        assert!(
            (actual_length - bond_length).abs() < 1e-9,
            "bond length should be exact by construction, got {actual_length}"
        );

        // Fully determined by the single real constraint (angle to CA),
        // exact regardless of the reference-driven azimuth choice.
        let angle_to_ca = angle_degrees(h, n, ca);
        assert!(
            (angle_to_ca - 120.0).abs() < 1e-6,
            "N-H should sit at the exact sp2 angle from CA by construction, got {angle_to_ca}"
        );

        // The actual claim this function makes over the arbitrary-azimuth
        // version: H's azimuth (viewed along the CA-N axis) should be on
        // the *opposite* side from the real previous-residue carbonyl
        // carbon, not the same side or perpendicular to it.
        let bond1 = normalize(subtract(ca, n));
        let perp_to_h = perpendicular_component(subtract(h, n), bond1);
        let perp_to_prev_c = perpendicular_component(subtract(prev_c, ca), bond1);
        let cos_azimuth = dot(perp_to_h, perp_to_prev_c);
        assert!(
            cos_azimuth < -0.9,
            "H's azimuth should be nearly opposite the previous residue's C, got cos = {cos_azimuth}"
        );
    }

    #[test]
    fn places_asparagines_amide_nh2_from_real_1ubq_geometry() {
        // Real ASN25 heavy-atom positions from 1UBQ. ND2 is bonded only
        // to CG; OD1 (the amide's other substituent on CG, part of the
        // same conjugated plane) fixes the orientation for HD21/HD22.
        let (pdb, _errors) = pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif")
            .expect("1UBQ fixture should parse");
        let asn25 = |name: &str| -> Point {
            pdb.atoms_with_hierarchy()
                .find(|h| {
                    h.residue().id().0 == 25
                        && h.residue().name() == Some("ASN")
                        && h.atom().name() == name
                })
                .unwrap_or_else(|| panic!("ASN25 {name} should be present"))
                .atom()
                .pos()
        };
        let nd2 = asn25("ND2");
        let cg = asn25("CG");
        let od1 = asn25("OD1");

        let bond_length = 1.02; // N-H: covalent radii 0.71 + 0.31
        let (hd21, hd22) = place_planar_amine_hydrogens(nd2, cg, od1, bond_length);

        for h in [hd21, hd22] {
            let actual_length = length(subtract(h, nd2));
            assert!(
                (actual_length - bond_length).abs() < 1e-9,
                "bond length should be exact by construction, got {actual_length}"
            );
            let angle_to_cg = angle_degrees(h, nd2, cg);
            assert!(
                (angle_to_cg - 120.0).abs() < 1e-6,
                "N-H should sit at the exact sp2 angle from CG by construction, got {angle_to_cg}"
            );
        }

        let angle_between_hydrogens = angle_degrees(hd21, nd2, hd22);
        assert!(
            (angle_between_hydrogens - 120.0).abs() < 1e-6,
            "the two hydrogens should also be 120 degrees apart from each other, got {angle_between_hydrogens}"
        );

        // The two hydrogens should be on opposite sides of the CG-ND2
        // axis (mirror images), not coincide or land on the same side.
        assert!(
            length(subtract(hd21, hd22)) > 1.0,
            "HD21 and HD22 should be clearly distinct positions"
        );

        // All four points (CG, ND2, HD21, HD22) should be coplanar --
        // that's the whole point of using OD1 to orient a shared plane.
        let bond1 = normalize(subtract(cg, nd2));
        let normal = {
            let a = subtract(hd21, nd2);
            let b = subtract(hd22, nd2);
            normalize((
                a.1 * b.2 - a.2 * b.1,
                a.2 * b.0 - a.0 * b.2,
                a.0 * b.1 - a.1 * b.0,
            ))
        };
        assert!(
            dot(bond1, normal).abs() < 1e-9,
            "CG should lie in the HD21-ND2-HD22 plane"
        );
    }

    #[test]
    fn places_serines_cb_methylene_hydrogens_from_real_1ubq_geometry() {
        // Real SER20 heavy-atom positions from 1UBQ. CB is bonded to CA
        // and OG (exactly two heavy neighbors) -- HB2/HB3 are the missing
        // methylene hydrogens this places.
        let (pdb, _errors) = pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif")
            .expect("1UBQ fixture should parse");
        let ser20 = |name: &str| -> Point {
            pdb.atoms_with_hierarchy()
                .find(|h| {
                    h.residue().id().0 == 20
                        && h.residue().name() == Some("SER")
                        && h.atom().name() == name
                })
                .unwrap_or_else(|| panic!("SER20 {name} should be present"))
                .atom()
                .pos()
        };
        let cb = ser20("CB");
        let ca = ser20("CA");
        let og = ser20("OG");

        let bond_length = 1.07; // C-H: covalent radii 0.76 + 0.31
        let (hb2, hb3) = place_tetrahedral_pair_from_two_neighbors(cb, ca, og, bond_length);

        for h in [hb2, hb3] {
            let actual_length = length(subtract(h, cb));
            assert!(
                (actual_length - bond_length).abs() < 1e-9,
                "bond length should be exact by construction, got {actual_length}"
            );
        }

        // Guaranteed exactly by construction, regardless of how close the
        // real CA-CB-OG angle is to ideal tetrahedral (it's ~110-111
        // degrees in practice, not exactly 109.47 -- see the CA/HA test's
        // note on real backbone geometry).
        let angle_between_hydrogens = angle_degrees(hb2, cb, hb3);
        assert!(
            (angle_between_hydrogens - 109.471_22).abs() < 1e-6,
            "HB2-CB-HB3 should be exactly tetrahedral by construction, got {angle_between_hydrogens}"
        );

        // Real, independent sanity bound on the two angles that aren't
        // exactly guaranteed (they depend on the real, slightly-distorted
        // CA-CB-OG input angle).
        for neighbor in [ca, og] {
            for h in [hb2, hb3] {
                let angle = angle_degrees(h, cb, neighbor);
                assert!(
                    (95.0..125.0).contains(&angle),
                    "HB-CB-neighbor angle {angle} is nowhere near a plausible tetrahedral range"
                );
            }
        }

        assert!(
            length(subtract(hb2, hb3)) > 1.0,
            "HB2 and HB3 should be clearly distinct positions"
        );
    }

    #[test]
    fn places_methionines_ce_methyl_hydrogens_from_real_1ubq_geometry() {
        // Real MET1 heavy-atom positions from 1UBQ. CE is bonded only to
        // SD -- HE1/HE2/HE3 are the missing terminal methyl hydrogens.
        let (pdb, _errors) = pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif")
            .expect("1UBQ fixture should parse");
        let met1 = |name: &str| -> Point {
            pdb.atoms_with_hierarchy()
                .find(|h| h.residue().id().0 == 1 && h.atom().name() == name)
                .unwrap_or_else(|| panic!("MET1 {name} should be present"))
                .atom()
                .pos()
        };
        let ce = met1("CE");
        let sd = met1("SD");

        let bond_length = 1.07; // C-H: covalent radii 0.76 + 0.31
        let (he1, he2, he3) = place_methyl_hydrogens(ce, sd, bond_length);
        let hydrogens = [he1, he2, he3];

        for h in hydrogens {
            let actual_length = length(subtract(h, ce));
            assert!(
                (actual_length - bond_length).abs() < 1e-9,
                "bond length should be exact by construction, got {actual_length}"
            );
            // Guaranteed exactly by construction: the single real
            // constraint (angle to SD) is fixed regardless of azimuth.
            let angle_to_sd = angle_degrees(h, ce, sd);
            assert!(
                (angle_to_sd - 109.471_22).abs() < 1e-6,
                "H-CE-SD should sit at the exact tetrahedral angle by construction, got {angle_to_sd}"
            );
        }

        // Also guaranteed exactly by construction for an ideal methyl: the
        // three H-CE-H angles equal the same tetrahedral angle too (a
        // real property of a symmetric tetrahedral arrangement, not just
        // a coincidence of this formula).
        for (a, b) in [(he1, he2), (he2, he3), (he1, he3)] {
            let angle = angle_degrees(a, ce, b);
            assert!(
                // Slightly looser than the direct angle-to-SD check above:
                // this one compounds more trig operations (two azimuth
                // rotations instead of one), so floating-point error
                // accumulates a bit further.
                (angle - 109.471_22).abs() < 1e-4,
                "H-CE-H should also be exactly tetrahedral by construction, got {angle}"
            );
        }

        // All three should be distinct (not collapsed onto each other by
        // an azimuth bug).
        assert!(length(subtract(he1, he2)) > 1.0);
        assert!(length(subtract(he2, he3)) > 1.0);
        assert!(length(subtract(he1, he3)) > 1.0);
    }
}
