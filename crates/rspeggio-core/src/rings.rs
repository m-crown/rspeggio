// crates/rspeggio-core/src/rings.rs
//
// M6: aromatic ring perception and plane geometry, built the same way
// `typing.rs` handles the rest of atom chemistry -- exact traversal over
// the CCD's own explicit aromatic bonds, no geometric/SMARTS re-perception
// (decisions 01/02). This is a two-stage split, same shape as
// `typing.rs`/`features.rs`:
//   - `perceive_rings`: pure chemistry-graph work on a `CcdComponent`, no
//     structure instance involved -- the same ring atom-id set for every
//     residue of that type, wherever it appears.
//   - `ring_geometry`: takes one perceived ring's atom ids plus a real
//     `Residue`'s resolved positions and computes the actual center/normal
//     for that specific instance.
//
// Ring perception: for each aromatic bond (u, v), the shortest path from u
// to v *not* using that bond, found by BFS over the aromatic-bond-only
// subgraph, is the smallest ring containing that bond -- closing the path
// with the excluded edge gives the cycle. Applying this per-edge and
// deduping by atom-id set recovers a correct minimal ring set even for a
// fused bicyclic system: tryptophan's indole has a shared CD2-CE2 edge
// between its 5- and 6-membered rings, and the edges unique to each ring
// only have a shortest path back through their own ring (confirmed against
// TRP's real bond graph in this module's tests).
//
// Ring geometry (center, normal) uses the same construction real
// pdbe-arpeggio gets from OpenBabel's `findCenterAndNormal` -- centroid of
// the ring atoms, and Newell's method for the normal (a sum of
// cross-product terms around the ring's traversal order, robust to a ring
// that isn't perfectly planar in real coordinates, unlike a single
// three-point cross product).
//
// `classify_ring_atom` covers the other half of real pdbe-arpeggio's ring
// contact code (`__calculate_atom_plane_contacts`): a non-aromatic atom
// sitting close to a ring's face -- cation-pi, donor-pi, carbon-pi
// (weak-donor CH...pi), and methionine-sulfur-pi. `HALOGENPI` isn't
// covered: it needs "xbond donor" atom typing, which doesn't exist yet
// (same gap the handoff notes for `XBOND` generally).

use crate::typing::AtomTypeBits;
use pdbtbx::Residue;
use rspeggio_ccd::component::CcdComponent;
use std::collections::{HashMap, HashSet, VecDeque};

type Point = (f64, f64, f64);

fn subtract(a: Point, b: Point) -> Point {
    (a.0 - b.0, a.1 - b.1, a.2 - b.2)
}

fn length(v: Point) -> f64 {
    (v.0 * v.0 + v.1 * v.1 + v.2 * v.2).sqrt()
}

fn normalize(v: Point) -> Point {
    let len = length(v);
    (v.0 / len, v.1 / len, v.2 / len)
}

fn dot(a: Point, b: Point) -> f64 {
    a.0 * b.0 + a.1 * b.1 + a.2 * b.2
}

// A perceived aromatic ring: just the chemistry-graph shape, in cyclic
// traversal order (needed for `ring_geometry`'s Newell's-method normal,
// which depends on going around the ring rather than an arbitrary order).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RingAtoms {
    pub atom_ids: Vec<String>,
}

// Builds an adjacency list over only the component's aromatic bonds.
fn aromatic_adjacency(component: &CcdComponent) -> HashMap<&str, Vec<&str>> {
    let mut adjacency: HashMap<&str, Vec<&str>> = HashMap::new();
    for bond in component.bonds() {
        if !bond.aromatic() {
            continue;
        }
        adjacency
            .entry(bond.atom_id_1())
            .or_default()
            .push(bond.atom_id_2());
        adjacency
            .entry(bond.atom_id_2())
            .or_default()
            .push(bond.atom_id_1());
    }
    adjacency
}

// Shortest path from `start` to `goal` over `adjacency`, never taking the
// single direct `start`-`goal` step (that's the bond being closed into a
// ring, not part of the path around it). `None` if `goal` is unreachable
// without it.
fn shortest_path_excluding_direct_edge<'a>(
    adjacency: &HashMap<&'a str, Vec<&'a str>>,
    start: &'a str,
    goal: &'a str,
) -> Option<Vec<&'a str>> {
    let mut visited: HashSet<&str> = HashSet::from([start]);
    let mut parent: HashMap<&str, &str> = HashMap::new();
    let mut queue: VecDeque<&str> = VecDeque::from([start]);

    while let Some(current) = queue.pop_front() {
        for &next in adjacency.get(current).into_iter().flatten() {
            if current == start && next == goal {
                continue; // the direct edge itself, not a path around it
            }
            if !visited.insert(next) {
                continue;
            }
            parent.insert(next, current);
            if next == goal {
                let mut path = vec![goal];
                let mut node = goal;
                while let Some(&p) = parent.get(node) {
                    path.push(p);
                    node = p;
                }
                path.reverse();
                return Some(path);
            }
            queue.push_back(next);
        }
    }
    None
}

// Perceives every aromatic ring in a component, as a minimal set of
// cycles over its aromatic-bond subgraph (see module doc). Pure chemistry,
// no structure instance -- the same result for every residue of this
// component's type.
pub fn perceive_rings(component: &CcdComponent) -> Vec<RingAtoms> {
    let adjacency = aromatic_adjacency(component);
    let mut seen: HashSet<Vec<String>> = HashSet::new();
    let mut rings = Vec::new();

    for bond in component.bonds() {
        if !bond.aromatic() {
            continue;
        }
        let Some(path) =
            shortest_path_excluding_direct_edge(&adjacency, bond.atom_id_1(), bond.atom_id_2())
        else {
            continue;
        };
        let mut key: Vec<String> = path.iter().map(|s| s.to_string()).collect();
        key.sort();
        if seen.insert(key) {
            rings.push(RingAtoms {
                atom_ids: path.into_iter().map(|s| s.to_string()).collect(),
            });
        }
    }
    rings
}

// One perceived ring's real geometry in a specific residue instance.
pub struct RingGeometry {
    pub center: Point,
    pub normal: Point, // unit vector; direction (vs `normal_opp`) is arbitrary
}

// Resolves `ring`'s real center and normal from `residue`'s actual atom
// positions. `None` if any ring atom isn't resolved in this residue (a
// structure that's missing a ring atom can't have its plane computed at
// all -- there's no partial-ring geometry to fall back to).
pub fn ring_geometry(ring: &RingAtoms, residue: &Residue) -> Option<RingGeometry> {
    let positions: Vec<Point> = ring
        .atom_ids
        .iter()
        .map(|id| residue.atoms().find(|a| a.name() == id).map(|a| a.pos()))
        .collect::<Option<Vec<_>>>()?;

    let n = positions.len() as f64;
    let sum = positions.iter().fold((0.0, 0.0, 0.0), |acc, p| {
        (acc.0 + p.0, acc.1 + p.1, acc.2 + p.2)
    });
    let center = (sum.0 / n, sum.1 / n, sum.2 / n);

    let mut normal = (0.0, 0.0, 0.0);
    for i in 0..positions.len() {
        let p1 = positions[i];
        let p2 = positions[(i + 1) % positions.len()];
        normal.0 += (p1.1 - p2.1) * (p1.2 + p2.2);
        normal.1 += (p1.2 - p2.2) * (p1.0 + p2.0);
        normal.2 += (p1.0 - p2.0) * (p1.1 + p2.1);
    }

    Some(RingGeometry {
        center,
        normal: normalize(normal),
    })
}

// The angle between two axes (ring normals, or a normal and a displacement
// vector), folded into [0, 90] degrees -- an axis has no inherent
// direction (a ring's `normal` vs `normal_opp` is an arbitrary choice), so
// the meaningful angle between two axes is always the acute one. Matches
// real pdbe-arpeggio's `group_angle`/`group_group_angle` with
// `signed=True` followed by `abs()` at every call site (`interactions.py`
// ring-ring/ring-atom code): mapping the raw arccos result into
// (-90, 90] and taking its absolute value is exactly this fold.
fn axis_angle_degrees(a: Point, b: Point) -> f64 {
    let cos_angle = dot(normalize(a), normalize(b)).clamp(-1.0, 1.0);
    let degrees = cos_angle.acos().to_degrees();
    if degrees > 90.0 {
        180.0 - degrees
    } else {
        degrees
    }
}

// The 9 ring-ring (plane-plane) geometric interaction types real
// pdbe-arpeggio classifies by binning `dihedral` (angle between the two
// ring planes) and `theta` (angle between one ring's plane and the
// center-to-center displacement) each into <=30/<=60/<=90 degree tiers --
// see `interactions.py`'s `__calculate_plane_plane_contacts`. Loosely:
// low dihedral + low theta is face-to-face stacking (`Ff`), low dihedral +
// high theta is offset/edge stacking (`Ee`), high dihedral + low theta is
// T-shaped/edge-to-face (`Fe`), and so on -- the naming itself (F/O/E for
// dihedral tier, f/t/e for theta tier) is real pdbe-arpeggio's own
// convention, kept here for the same reason CCD names are kept verbatim:
// so this stays checkable against the oracle.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RingRingInteraction {
    Ff,
    Of,
    Ee,
    Ft,
    Ot,
    Et,
    Fe,
    Oe,
    Ef,
}

// Real pdbe-arpeggio's ring-ring centroid-distance cutoff before even
// computing the angles (`config.CONTACT_TYPES['aromatic']['centroid_distance']`).
pub const CENTROID_DISTANCE_MAX: f64 = 6.0;

// Classifies the geometric relationship between two real ring instances.
// `None` if they're farther apart than `CENTROID_DISTANCE_MAX` -- not
// itself an aromatic-stacking judgement (that's a caller's job, same as
// `contacts::classify_distance` not deciding which pairs are chemically
// interesting), just the raw geometry.
pub fn classify_ring_ring(
    a: &RingGeometry,
    b: &RingGeometry,
) -> Option<(f64, RingRingInteraction)> {
    let distance = length(subtract(a.center, b.center));
    if distance > CENTROID_DISTANCE_MAX {
        return None;
    }

    let dihedral = axis_angle_degrees(a.normal, b.normal);
    let theta = axis_angle_degrees(a.normal, subtract(a.center, b.center));

    // Both `dihedral` and `theta` are already folded into [0, 90] by
    // `axis_angle_degrees`, so each tier's final branch (theta/dihedral
    // <= 90) is exhaustive -- not a fallback for an out-of-range value.
    // Deliberately *not* refactored into an orthogonal per-axis formula:
    // the real 9-label naming isn't a clean product of independent
    // per-tier letters (the `theta<=90` row reads EE/ET/EF, the reverse
    // letter order of the `theta<=30`/`theta<=60` rows' F/T/E), so mirroring
    // real pdbe-arpeggio's own if/elif chain verbatim (`interactions.py`'s
    // `__calculate_plane_plane_contacts`) is both simpler and less error
    // prone than reverse-engineering a formula for it.
    let kind = if dihedral <= 30.0 {
        if theta <= 30.0 {
            RingRingInteraction::Ff
        } else if theta <= 60.0 {
            RingRingInteraction::Of
        } else {
            RingRingInteraction::Ee
        }
    } else if dihedral <= 60.0 {
        if theta <= 30.0 {
            RingRingInteraction::Ft
        } else if theta <= 60.0 {
            RingRingInteraction::Ot
        } else {
            RingRingInteraction::Et
        }
    } else if theta <= 30.0 {
        RingRingInteraction::Fe
    } else if theta <= 60.0 {
        RingRingInteraction::Oe
    } else {
        RingRingInteraction::Ef
    };

    Some((distance, kind))
}

// The 4 (of 5 real) ring-atom contact types this project can currently
// type. Real pdbe-arpeggio also has `HALOGENPI` (a halogen-bond-donor
// atom near a ring face) -- omitted here, not silently: it needs "xbond
// donor" atom typing this project hasn't built (see module doc).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RingAtomInteraction {
    CarbonPi,
    CationPi,
    DonorPi,
    MetSulphurPi,
}

// Real pdbe-arpeggio's ring-atom thresholds
// (`config.CONTACT_TYPES['aromatic']['atom_aromatic_distance']`/
// `['met_sulphur_aromatic_distance']`). `MET_SULPHUR_DISTANCE_MAX` is
// deliberately much larger (6.0A vs 4.5A) -- real pdbe-arpeggio gives a
// methionine sulfur's lone pair a longer reach toward a ring face than an
// ordinary donor/cation, reflecting sulfur's larger, more polarizable
// electron cloud.
pub const ATOM_AROMATIC_DISTANCE_MAX: f64 = 4.5;
pub const MET_SULPHUR_DISTANCE_MAX: f64 = 6.0;

// Classifies every ring-atom interaction between `ring` and one nearby
// atom, given that atom's real position, element, owning residue name, and
// already-computed `AtomTypeBits` (from `typing::type_atom` -- this
// function takes bits rather than a `CcdAtom`+`CcdComponent` pair so it
// stays agnostic to how the caller got them, same as `features.rs`'s
// `classify_features` takes already-typed atoms rather than re-deriving
// types itself). Can return more than one interaction (e.g. a weak-donor
// carbon that's also positively ionisable would be both CARBONPI and
// CATIONPI, mirroring real pdbe-arpeggio's `potential_interactions` set).
pub fn classify_ring_atom(
    ring: &RingGeometry,
    atom_pos: Point,
    atom_element: &str,
    atom_residue_name: &str,
    atom_bits: AtomTypeBits,
) -> Vec<RingAtomInteraction> {
    let mut interactions = Vec::new();

    // No aromatic-atom-to-ring interactions -- an aromatic atom belongs to
    // its own ring's plane-plane classification instead (`classify_ring_ring`).
    if atom_bits.contains(AtomTypeBits::AROMATIC) {
        return interactions;
    }

    let distance = length(subtract(atom_pos, ring.center));

    if distance <= ATOM_AROMATIC_DISTANCE_MAX {
        let theta = axis_angle_degrees(ring.normal, subtract(ring.center, atom_pos));
        if theta <= 30.0 {
            if atom_element.eq_ignore_ascii_case("C")
                && atom_bits.contains(AtomTypeBits::WEAK_HBOND_DONOR)
            {
                interactions.push(RingAtomInteraction::CarbonPi);
            }
            if atom_bits.contains(AtomTypeBits::POS_IONISABLE) {
                interactions.push(RingAtomInteraction::CationPi);
            }
            if atom_bits.contains(AtomTypeBits::HBOND_DONOR) {
                interactions.push(RingAtomInteraction::DonorPi);
            }
        }
    }

    if distance <= MET_SULPHUR_DISTANCE_MAX
        && atom_element.eq_ignore_ascii_case("S")
        && atom_residue_name.eq_ignore_ascii_case("MET")
    {
        interactions.push(RingAtomInteraction::MetSulphurPi);
    }

    interactions
}

#[cfg(test)]
mod tests {
    use super::*;
    use rspeggio_ccd::common::common_components;

    fn ring_atom_set(component: &CcdComponent) -> Vec<HashSet<String>> {
        perceive_rings(component)
            .into_iter()
            .map(|r| r.atom_ids.into_iter().collect())
            .collect()
    }

    #[test]
    fn phenylalanine_has_exactly_one_six_membered_ring() {
        let components = common_components();
        let phe = components.get("PHE").expect("PHE should be bundled");

        let rings = ring_atom_set(phe);
        assert_eq!(rings.len(), 1, "PHE's side chain is a single benzene ring");

        let expected: HashSet<String> = ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        assert_eq!(rings[0], expected);
    }

    #[test]
    fn histidine_has_exactly_one_five_membered_ring() {
        let components = common_components();
        let his = components.get("HIS").expect("HIS should be bundled");

        let rings = ring_atom_set(his);
        assert_eq!(
            rings.len(),
            1,
            "HIS's side chain is a single imidazole ring"
        );

        let expected: HashSet<String> = ["CG", "ND1", "CD2", "CE1", "NE2"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        assert_eq!(rings[0], expected);
    }

    #[test]
    fn tryptophan_has_two_fused_rings() {
        // Indole: a 5-membered pyrrole ring fused to a 6-membered benzo
        // ring, sharing the CD2-CE2 edge -- the real case that needs
        // per-edge shortest-cycle perception rather than "one connected
        // aromatic component = one ring".
        let components = common_components();
        let trp = components.get("TRP").expect("TRP should be bundled");

        let rings = ring_atom_set(trp);
        assert_eq!(rings.len(), 2, "TRP's indole has two fused rings");

        let five: HashSet<String> = ["CG", "CD1", "NE1", "CE2", "CD2"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let six: HashSet<String> = ["CD2", "CE2", "CZ2", "CH2", "CZ3", "CE3"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        assert!(
            rings.contains(&five),
            "expected the 5-membered pyrrole ring, got {rings:?}"
        );
        assert!(
            rings.contains(&six),
            "expected the 6-membered benzo ring, got {rings:?}"
        );
    }

    #[test]
    fn alanine_has_no_rings() {
        let components = common_components();
        let ala = components.get("ALA").expect("ALA should be bundled");
        assert!(perceive_rings(ala).is_empty());
    }

    #[test]
    fn a_real_phenylalanine_rings_geometry_has_a_normal_roughly_perpendicular_to_every_ring_bond() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");
        let components = common_components();
        let phe = components.get("PHE").expect("PHE should be bundled");
        let ring = perceive_rings(phe)
            .into_iter()
            .next()
            .expect("PHE has a ring");

        let residue = pdb
            .residues()
            .find(|r| r.name() == Some("PHE"))
            .expect("1UBQ has at least one PHE");

        let geometry = ring_geometry(&ring, residue).expect("all ring atoms should be resolved");

        // A real (near-)planar aromatic ring: every bond vector within the
        // ring should lie close to the plane, i.e. be close to
        // perpendicular to the computed normal.
        let positions: Vec<Point> = ring
            .atom_ids
            .iter()
            .map(|id| residue.atoms().find(|a| a.name() == id).unwrap().pos())
            .collect();
        for i in 0..positions.len() {
            let next = positions[(i + 1) % positions.len()];
            let bond_vector = normalize(subtract(next, positions[i]));
            let cos_to_normal = dot(bond_vector, geometry.normal).abs();
            assert!(
                cos_to_normal < 0.2,
                "ring bond {i} should be close to perpendicular to the normal, got cos = {cos_to_normal}"
            );
        }

        // The center should sit close to the average position, well within
        // the ring's own extent from any single atom.
        for p in &positions {
            assert!(length(subtract(*p, geometry.center)) < 2.0);
        }
    }

    #[test]
    fn two_coplanar_rings_far_apart_are_face_to_face() {
        // Two hand-built hexagonal rings stacked directly on top of each
        // other, same orientation, 3.5A apart -- the textbook Ff
        // (face-to-face) pi-stack geometry: dihedral 0, theta 0.
        let hexagon = |z: f64| -> Vec<Point> {
            (0..6)
                .map(|i| {
                    let angle = std::f64::consts::PI / 3.0 * i as f64;
                    (angle.cos() * 1.4, angle.sin() * 1.4, z)
                })
                .collect()
        };
        // Build geometry directly rather than through a real Residue --
        // this is pure hand-built Euclidean geometry, not chemistry, so a
        // synthetic case is the right (only) way to pin an exact expected
        // classification down.
        let geometry_from_points = |points: Vec<Point>| -> RingGeometry {
            let n = points.len() as f64;
            let sum = points.iter().fold((0.0, 0.0, 0.0), |acc, p| {
                (acc.0 + p.0, acc.1 + p.1, acc.2 + p.2)
            });
            let center = (sum.0 / n, sum.1 / n, sum.2 / n);
            let mut normal = (0.0, 0.0, 0.0);
            for i in 0..points.len() {
                let p1 = points[i];
                let p2 = points[(i + 1) % points.len()];
                normal.0 += (p1.1 - p2.1) * (p1.2 + p2.2);
                normal.1 += (p1.2 - p2.2) * (p1.0 + p2.0);
                normal.2 += (p1.0 - p2.0) * (p1.1 + p2.1);
            }
            RingGeometry {
                center,
                normal: normalize(normal),
            }
        };

        let a = geometry_from_points(hexagon(0.0));
        let b = geometry_from_points(hexagon(3.5));

        let (distance, kind) = classify_ring_ring(&a, &b).expect("well within centroid distance");
        assert!((distance - 3.5).abs() < 1e-9);
        assert_eq!(kind, RingRingInteraction::Ff);
    }

    #[test]
    fn two_perpendicular_rings_are_edge_to_face() {
        // Ring B rotated 90 degrees relative to A, offset along A's own
        // plane so its center is roughly in-plane with A -- the textbook
        // T-shaped/edge-to-face geometry: dihedral ~90, theta ~90.
        let hexagon_xy = |center: Point| -> Vec<Point> {
            (0..6)
                .map(|i| {
                    let angle = std::f64::consts::PI / 3.0 * i as f64;
                    (
                        center.0 + angle.cos() * 1.4,
                        center.1 + angle.sin() * 1.4,
                        center.2,
                    )
                })
                .collect()
        };
        let hexagon_xz = |center: Point| -> Vec<Point> {
            (0..6)
                .map(|i| {
                    let angle = std::f64::consts::PI / 3.0 * i as f64;
                    (
                        center.0 + angle.cos() * 1.4,
                        center.1,
                        center.2 + angle.sin() * 1.4,
                    )
                })
                .collect()
        };
        let geometry_from_points = |points: Vec<Point>| -> RingGeometry {
            let n = points.len() as f64;
            let sum = points.iter().fold((0.0, 0.0, 0.0), |acc, p| {
                (acc.0 + p.0, acc.1 + p.1, acc.2 + p.2)
            });
            let center = (sum.0 / n, sum.1 / n, sum.2 / n);
            let mut normal = (0.0, 0.0, 0.0);
            for i in 0..points.len() {
                let p1 = points[i];
                let p2 = points[(i + 1) % points.len()];
                normal.0 += (p1.1 - p2.1) * (p1.2 + p2.2);
                normal.1 += (p1.2 - p2.2) * (p1.0 + p2.0);
                normal.2 += (p1.0 - p2.0) * (p1.1 + p2.1);
            }
            RingGeometry {
                center,
                normal: normalize(normal),
            }
        };

        let a = geometry_from_points(hexagon_xy((0.0, 0.0, 0.0)));
        // B's plane contains the Z axis (A's normal), and B is offset
        // along +X from A within A's own plane -- so the displacement
        // A->B is perpendicular to A's normal (theta ~90), and B's normal
        // (Y axis) is perpendicular to A's normal (Z axis) too (dihedral
        // ~90).
        let b = geometry_from_points(hexagon_xz((4.5, 0.0, 0.0)));

        let (_, kind) = classify_ring_ring(&a, &b).expect("well within centroid distance");
        assert_eq!(kind, RingRingInteraction::Ef);
    }

    #[test]
    fn a_real_ring_ring_contact_is_found_between_two_phenylalanines_in_1ca2() {
        // Found by an ad-hoc sweep of every ring pair in 1CA2 (which has 45
        // real perceived ring instances): PHE66 and PHE95 sit ~5.4A apart
        // centroid-to-centroid, well within the 6A cutoff -- confirms
        // `ring_geometry` + `classify_ring_ring` work end to end on real,
        // independently-resolved structure coordinates, not just hand-built
        // synthetic hexagons.
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1CA2.cif").expect("1CA2 should load");
        let components = common_components();
        let phe = components.get("PHE").expect("PHE should be bundled");
        let ring = perceive_rings(phe)
            .into_iter()
            .next()
            .expect("PHE has a ring");

        let residue_named = |seq_id: isize| {
            pdb.residues()
                .find(|r| r.id().0 == seq_id && r.name() == Some("PHE"))
                .unwrap_or_else(|| panic!("PHE{seq_id} should be present in 1CA2"))
        };
        let phe66 = ring_geometry(&ring, residue_named(66)).expect("PHE66's ring should resolve");
        let phe95 = ring_geometry(&ring, residue_named(95)).expect("PHE95's ring should resolve");

        let (distance, _kind) =
            classify_ring_ring(&phe66, &phe95).expect("PHE66/PHE95 should be within 6A");
        assert!(
            (5.0..6.0).contains(&distance),
            "expected PHE66-PHE95 centroid distance around 5.4A, got {distance}"
        );
    }

    // Shared by every `classify_ring_atom` real-data test below: BPTI
    // (5PTI) is small and already a fixture elsewhere in this crate, and
    // an ad-hoc sweep of every ring-atom pair in it (plus 1UBQ/1MBO/1FLV/
    // 4FXC/1CA2) turned up real, independently-verifiable examples of
    // every interaction this project can currently type.
    fn bpti_ring(comp_id: &str, components: &HashMap<String, CcdComponent>) -> RingAtoms {
        let component = components
            .get(comp_id)
            .unwrap_or_else(|| panic!("{comp_id} should be bundled"));
        perceive_rings(component)
            .into_iter()
            .next()
            .unwrap_or_else(|| panic!("{comp_id} should have a ring"))
    }

    fn bpti_residue<'a>(pdb: &'a pdbtbx::PDB, seq_id: isize, comp_id: &str) -> &'a pdbtbx::Residue {
        pdb.residues()
            .find(|r| r.id().0 == seq_id && r.name() == Some(comp_id))
            .unwrap_or_else(|| panic!("{comp_id}{seq_id} should be present in BPTI"))
    }

    #[test]
    fn a_real_carbon_pi_contact_is_found_in_bpti() {
        // PHE4's ring face sits close to ARG42's CB -- a plain aliphatic
        // carbon with an attached hydrogen (weak hbond donor), not
        // ionisable or a strong donor, so CARBONPI should be the only
        // interaction found.
        let (pdb, _errors) =
            pdbtbx::open("tests/fixtures/structures/5PTI.cif").expect("BPTI should load");
        let components = common_components();

        let ring = bpti_ring("PHE", &components);
        let phe4 = bpti_residue(&pdb, 4, "PHE");
        let geometry = ring_geometry(&ring, phe4).expect("PHE4's ring should resolve");

        let arg42 = bpti_residue(&pdb, 42, "ARG");
        let cb = arg42
            .atoms()
            .find(|a| a.name() == "CB")
            .expect("ARG42 should have a CB");
        let arg_component = components.get("ARG").expect("ARG should be bundled");
        let cb_ccd = arg_component
            .atoms()
            .iter()
            .find(|a| a.atom_id() == "CB")
            .expect("ARG's CCD entry should have a CB");
        let bits = crate::typing::type_atom(cb_ccd, arg_component);

        let interactions = classify_ring_atom(&geometry, cb.pos(), "C", "ARG", bits);
        assert_eq!(
            interactions,
            vec![RingAtomInteraction::CarbonPi],
            "expected only CARBONPI between PHE4's ring and ARG42 CB, got {interactions:?}"
        );
    }

    #[test]
    fn a_real_cation_pi_and_donor_pi_contact_is_found_in_bpti() {
        // TYR10's ring sits close to LYS41's real, fully protonated NZ
        // (NH3+) -- positively ionisable *and* a real donor, so both
        // CATIONPI and DONORPI should fire for the same atom.
        let (pdb, _errors) =
            pdbtbx::open("tests/fixtures/structures/5PTI.cif").expect("BPTI should load");
        let components = common_components();

        let ring = bpti_ring("TYR", &components);
        let tyr10 = bpti_residue(&pdb, 10, "TYR");
        let geometry = ring_geometry(&ring, tyr10).expect("TYR10's ring should resolve");

        let lys41 = bpti_residue(&pdb, 41, "LYS");
        let nz = lys41
            .atoms()
            .find(|a| a.name() == "NZ")
            .expect("LYS41 should have an NZ");
        let lys_component = components.get("LYS").expect("LYS should be bundled");
        let nz_ccd = lys_component
            .atoms()
            .iter()
            .find(|a| a.atom_id() == "NZ")
            .expect("LYS's CCD entry should have an NZ");
        let bits = crate::typing::type_atom(nz_ccd, lys_component);

        let interactions = classify_ring_atom(&geometry, nz.pos(), "N", "LYS", bits);
        assert!(interactions.contains(&RingAtomInteraction::CationPi));
        assert!(interactions.contains(&RingAtomInteraction::DonorPi));
    }

    #[test]
    fn a_real_met_sulphur_pi_contact_is_found_in_bpti() {
        // TYR23's ring sits within the wider 6A sulfur-specific cutoff of
        // MET52's real SD -- the one ring-atom interaction that isn't
        // gated by the tighter 4.5A/30-degree "near the face" check at
        // all.
        let (pdb, _errors) =
            pdbtbx::open("tests/fixtures/structures/5PTI.cif").expect("BPTI should load");
        let components = common_components();

        let ring = bpti_ring("TYR", &components);
        let tyr23 = bpti_residue(&pdb, 23, "TYR");
        let geometry = ring_geometry(&ring, tyr23).expect("TYR23's ring should resolve");

        let met52 = bpti_residue(&pdb, 52, "MET");
        let sd = met52
            .atoms()
            .find(|a| a.name() == "SD")
            .expect("MET52 should have an SD");
        let met_component = components.get("MET").expect("MET should be bundled");
        let sd_ccd = met_component
            .atoms()
            .iter()
            .find(|a| a.atom_id() == "SD")
            .expect("MET's CCD entry should have an SD");
        let bits = crate::typing::type_atom(sd_ccd, met_component);

        let interactions = classify_ring_atom(&geometry, sd.pos(), "S", "MET", bits);
        assert_eq!(interactions, vec![RingAtomInteraction::MetSulphurPi]);
    }

    #[test]
    fn an_aromatic_atom_never_gets_a_ring_atom_interaction() {
        // The explicit "no aromatic atom-ring interactions" guard: even an
        // atom that would otherwise satisfy e.g. DONORPI/CARBONPI
        // shouldn't fire if it's itself typed aromatic (it belongs to its
        // own ring's plane-plane classification instead).
        let geometry = RingGeometry {
            center: (0.0, 0.0, 0.0),
            normal: (0.0, 0.0, 1.0),
        };

        let interactions = classify_ring_atom(
            &geometry,
            (1.0, 0.0, 0.0),
            "C",
            "PHE",
            AtomTypeBits::AROMATIC | AtomTypeBits::WEAK_HBOND_DONOR,
        );
        assert!(interactions.is_empty());
    }
}
