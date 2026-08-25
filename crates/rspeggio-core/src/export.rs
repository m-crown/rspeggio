// crates/rspeggio-core/src/export.rs
//
// M7: JSON export matching real pdbe-arpeggio's own schema
// (`Arpeggio.get_contacts` in `interactions.py`) -- atom-atom, plus the
// four ring/amide "plane" contact types `rings.rs` already computes the
// geometry for: plane-plane (ring-ring), atom-plane (ring-atom),
// group-group (amide-amide), group-plane (amide-ring).
//
// Real, confirmed limitation worth knowing before reading further: ring
// perception (`rings::perceive_rings`) only sees rings whose bonds are
// flagged `pdbx_aromatic_flag = Y` in their real CCD entry (decisions
// 01/04 -- no geometric re-perception). HEM's real CCD entry flags this
// *inconsistently* across its own four chemically-equivalent pyrrole
// rings (checked directly): rings A and C (`NA-C1A-C2A-C3A-C4A`,
// `NC-C1C-C2C-C3C-C4C`) are flagged aromatic and so are perceived here;
// rings B and D are not, and aren't. Real pdbe-arpeggio doesn't have this
// gap at all -- its OpenBabel-based aromaticity *re*-perception treats all
// four consistently (confirmed: `tests/fixtures/golden/1MBO.json` has a
// real CARBONPI atom-plane entry against ring B specifically, which this
// project's CCD-flag-based perception can't see). A real, principled
// consequence of the project's own no-re-perception decision, not a bug
// in the ring-perception algorithm itself -- but real CCD curation isn't
// perfectly self-consistent either, worth knowing precisely rather than
// assuming "HEM has no rings at all".
//
// Every export function here takes a `&selection::SelectionContext`
// (`selection.rs`): real pdbe-arpeggio's `-s` selection mechanism, which
// `interacting_entities` (INTER/INTRA_SELECTION/INTRA_NON_SELECTION/
// SELECTION_WATER/NON_SELECTION_WATER/WATER_WATER/INTRA_BINDING_SITE)
// depends on entirely -- and which real arpeggio also uses to restrict
// which contacts get computed at all (only pairs where both real atoms
// fall in the selection's binding-site expansion, `selection_plus`, are
// ever considered; see `selection.rs`'s own doc). `SelectionContext::
// whole_structure` treats every residue as selected, reproducing real
// arpeggio's own default behavior when `-s` is omitted -- 1UBQ's golden
// fixture really was generated this way (confirmed, not assumed: see
// `selection.rs`'s doc comment).
//
// This means most of this module's own real-data-verified tests DO now
// byte-match `tests/fixtures/golden/*.json`'s own `interacting_entities`
// values directly, using a real `RESNAME:<id>` selection worked out by
// hand-checking each golden fixture (not guessed) -- see the tests
// tagged "reproduces ... golden interacting_entities exactly" at the
// bottom of this module.

use crate::config::{self, DistanceCategory, FeatureBits};
use crate::contacts::{euclidean_distance, find_contacts, Contact};
use crate::features::classify_features;
use crate::rings::{
    self, AmideGeometry, AmideGroup, RingAtomInteraction, RingAtoms, RingGeometry,
    RingRingInteraction,
};
use pdbtbx::{
    AtomConformerResidueChainModel, ContainsAtomConformer, ContainsAtomConformerResidue,
    ContainsAtomConformerResidueChain, Residue, PDB,
};
use rayon::prelude::*;
use rspeggio_ccd::component::{CcdComponent, ComponentType};
use serde::Serialize;
use std::collections::HashMap;

#[derive(Debug, Clone, Serialize)]
pub struct AtomIdentity {
    pub auth_asym_id: String,
    pub auth_atom_id: String,
    pub auth_seq_id: isize,
    pub label_comp_id: String,
    pub label_comp_type: &'static str,
    #[serde(rename = "pdbx_PDB_ins_code")]
    pub pdbx_pdb_ins_code: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct AtomAtomContactJson {
    pub bgn: AtomIdentity,
    pub end: AtomIdentity,
    #[serde(rename = "type")]
    pub entry_type: &'static str,
    pub distance: f64,
    pub contact: Vec<&'static str>,
    pub interacting_entities: &'static str,
}

// The shape real `_prepare_plane_plane_contact_for_export` builds --
// shared by plane-plane (ring-ring), group-group (amide-amide), and
// group-plane (amide-ring), same as the oracle reuses one function for
// all three. `bgn`/`end` are `AtomIdentity`s built by `group_identity`
// (a ring/amide's own residue position, with `auth_atom_id` set to its
// member atoms rather than a single real atom).
#[derive(Debug, Clone, Serialize)]
pub struct PlanePlaneContactJson {
    pub bgn: AtomIdentity,
    pub end: AtomIdentity,
    #[serde(rename = "type")]
    pub entry_type: &'static str,
    pub distance: f64,
    pub contact: Vec<&'static str>,
    pub interacting_entities: &'static str,
}

// Real `_prepare_atom_plane_contact_for_export`'s shape (atom-plane,
// ring-atom): `bgn` is a real atom, `end` is the ring it's near.
#[derive(Debug, Clone, Serialize)]
pub struct AtomPlaneContactJson {
    pub bgn: AtomIdentity,
    pub end: AtomIdentity,
    #[serde(rename = "type")]
    pub entry_type: &'static str,
    pub distance: f64,
    pub contact: Vec<&'static str>,
    pub interacting_entities: &'static str,
}

// Real pdbe-arpeggio's own 15-position label list
// (`Arpeggio.get_contacts`, `interactions.py:178`) -- the first 5 are
// `DistanceCategory`, mutually exclusive (exactly one always present); the
// remaining 10 are `FeatureBits`, zero or more present. The feature order
// here matches `FeatureBits`' own declared bit order exactly (both are the
// same real SIFt position order), so this list doubles as a decode table
// for the bitflags themselves.
pub(crate) const FEATURE_LABELS: &[(FeatureBits, &str)] = &[
    (FeatureBits::HBOND, "hbond"),
    (FeatureBits::WEAK_HBOND, "weak_hbond"),
    (FeatureBits::XBOND, "xbond"),
    (FeatureBits::IONIC, "ionic"),
    (FeatureBits::METAL_COMPLEX, "metal_complex"),
    (FeatureBits::AROMATIC, "aromatic"),
    (FeatureBits::HYDROPHOBIC, "hydrophobic"),
    (FeatureBits::CARBONYL, "carbonyl"),
    (FeatureBits::POLAR, "polar"),
    (FeatureBits::WEAK_POLAR, "weak_polar"),
];

fn distance_category_label(category: DistanceCategory) -> &'static str {
    match category {
        DistanceCategory::Clash => "clash",
        DistanceCategory::Covalent => "covalent",
        DistanceCategory::VdwClash => "vdw_clash",
        DistanceCategory::Vdw => "vdw",
        DistanceCategory::Proximal => "proximal",
    }
}

fn contact_labels(category: DistanceCategory, features: FeatureBits) -> Vec<&'static str> {
    let mut labels = vec![distance_category_label(category)];
    for (bit, label) in FEATURE_LABELS {
        if features.contains(*bit) {
            labels.push(label);
        }
    }
    labels
}

// Real pdbe-arpeggio rounds with `round(np.float64(distance), 2)`
// (`interactions.py:190`), Python/numpy's round-half-to-even. This uses
// ordinary round-half-away-from-zero instead (`f64::round`) -- the two
// only disagree on an exact tie at the second decimal place, vanishingly
// unlikely for a real, continuously-varying Euclidean distance, but a real
// (if practically unobservable) divergence worth naming rather than
// silently assuming away.
fn round_2(x: f64) -> f64 {
    (x * 100.0).round() / 100.0
}

fn insertion_code(residue: &pdbtbx::Residue) -> String {
    match residue.id().1 {
        Some(code) => code.to_string(),
        None => " ".to_string(),
    }
}

fn atom_identity(atom: &AtomConformerResidueChainModel, component: &CcdComponent) -> AtomIdentity {
    let residue = atom.residue();
    AtomIdentity {
        auth_asym_id: atom.chain().id().to_string(),
        auth_atom_id: atom.atom().name().to_string(),
        auth_seq_id: residue.id().0,
        label_comp_id: residue.name().unwrap_or_default().to_string(),
        label_comp_type: component.component_type().code(),
        pdbx_pdb_ins_code: insertion_code(residue),
    }
}

// Real `__get_contact_type`'s logic verbatim (`interactions.py:643`):
// a sequence of unconditional overwrites, later ones winning, not an
// if/elif chain -- e.g. two real water atoms both pass the "both
// selected" check too when the selection is the whole structure, but
// `WATER_WATER` still wins because it's checked last. With `selection` a
// whole-structure `Selection` (`a_sel`/`b_sel` always true), this
// collapses to exactly 3 reachable outcomes (WATER_WATER/SELECTION_WATER/
// INTRA_SELECTION) -- confirmed against 1UBQ's own golden fixture, which
// really was generated with no selection at all.
fn interacting_entities_atom_atom(
    selection: &crate::selection::Selection,
    atom_a: &pdbtbx::Atom,
    component_a: ComponentType,
    atom_b: &pdbtbx::Atom,
    component_b: ComponentType,
) -> &'static str {
    let a_selected = selection.contains(atom_a);
    let b_selected = selection.contains(atom_b);
    let a_water = component_a == ComponentType::Water;
    let b_water = component_b == ComponentType::Water;

    let mut result = "";
    if !a_selected && !b_selected {
        result = "INTRA_NON_SELECTION";
    }
    if a_selected && b_selected {
        result = "INTRA_SELECTION";
    }
    if a_selected != b_selected {
        result = "INTER";
    }
    if (a_selected && b_water) || (b_selected && a_water) {
        result = "SELECTION_WATER";
    }
    if (!a_selected && b_water) || (!b_selected && a_water) {
        result = "NON_SELECTION_WATER";
    }
    if a_water && b_water {
        result = "WATER_WATER";
    }
    result
}

// Real ring-ring/atom-plane/group-group/group-plane contact-type logic
// (each `__calculate_*_contacts` function has its own copy, always this
// same shape): `strict_*` is real membership in the raw selection,
// `plus_*` in the binding-site-expanded one. Every caller here has
// already filtered to pairs where both sides pass `plus`, so
// `INTRA_NON_SELECTION` is unreachable in practice (kept anyway, to stay
// a faithful port of the real sequential-overwrite logic rather than a
// hand-simplified version of it).
fn interacting_entities_group(
    strict_a: bool,
    strict_b: bool,
    plus_a: bool,
    plus_b: bool,
) -> &'static str {
    let mut result = "";
    if !strict_a && !strict_b {
        result = "INTRA_NON_SELECTION";
    }
    if plus_a && plus_b {
        result = "INTRA_BINDING_SITE";
    }
    if strict_a && strict_b {
        result = "INTRA_SELECTION";
    }
    if strict_a != strict_b {
        result = "INTER";
    }
    result
}

fn component_for<'a>(
    atom: &AtomConformerResidueChainModel,
    components: &'a HashMap<String, CcdComponent>,
) -> Result<&'a CcdComponent, String> {
    let comp_id = atom.residue().name().unwrap_or_default();
    components
        .get(comp_id)
        .ok_or_else(|| format!("unknown CCD component: {comp_id}"))
}

fn build_atom_atom_entry(
    contact: &Contact,
    components: &HashMap<String, CcdComponent>,
    selection: &crate::selection::SelectionContext,
) -> Result<AtomAtomContactJson, String> {
    let component_1 = component_for(&contact.atom_1, components)?;
    let component_2 = component_for(&contact.atom_2, components)?;

    let features = classify_features(contact, components);
    let interacting_entities = interacting_entities_atom_atom(
        &selection.selection,
        contact.atom_1.atom(),
        component_1.component_type(),
        contact.atom_2.atom(),
        component_2.component_type(),
    );

    Ok(AtomAtomContactJson {
        bgn: atom_identity(&contact.atom_1, component_1),
        end: atom_identity(&contact.atom_2, component_2),
        entry_type: "atom-atom",
        distance: round_2(contact.distance),
        contact: contact_labels(contact.category, features),
        interacting_entities,
    })
}

// Every atom-atom contact in `pdb`, in real pdbe-arpeggio's JSON shape.
// Real arpeggio only ever computes contacts among `selection_plus` atoms
// at all (`_make_selection` rebuilds its neighbor search over exactly that
// set) -- both real atoms of a contact must pass that filter here too,
// not just be present in `pdb`. `Err` as soon as any *surviving* pair's
// residue comp_id has no matching `CcdComponent` -- consistent with
// decision 03 (unknown components fail loudly), since `label_comp_type`
// genuinely can't be produced without one.
pub fn export_atom_atom_contacts(
    pdb: &PDB,
    components: &HashMap<String, CcdComponent>,
    selection: &crate::selection::SelectionContext,
) -> Result<Vec<AtomAtomContactJson>, String> {
    find_contacts(pdb, config::CONTACT_TYPES_MAX_DIST)
        .par_iter()
        .filter(|c| {
            selection.selection_plus.contains(c.atom_1.atom())
                && selection.selection_plus.contains(c.atom_2.atom())
        })
        .map(|c| build_atom_atom_entry(c, components, selection))
        .collect()
}

// A "group" identity (a ring or amide's real position in the structure,
// e.g. `end.auth_atom_id` in real pdbe-arpeggio's own atom-plane JSON,
// `"C1B,C2B,C3B,C4B,NB"`) is the same shape as a single atom's
// (`AtomIdentity`) with `auth_atom_id` set to the sorted, comma-joined
// member atom names instead of one real atom's name -- real
// `_prepare_plane_plane_contact_for_export` builds it exactly this way
// (`make_pymol_json` on the *residue*, then separately overwrites
// `auth_atom_id`), so this reuses `AtomIdentity` rather than a separate type.
fn group_identity(
    hierarchy: &AtomConformerResidueChainModel,
    component: &CcdComponent,
    member_atom_ids: &[String],
) -> AtomIdentity {
    let residue = hierarchy.residue();
    let mut sorted_ids = member_atom_ids.to_vec();
    sorted_ids.sort();
    AtomIdentity {
        auth_asym_id: hierarchy.chain().id().to_string(),
        auth_atom_id: sorted_ids.join(","),
        auth_seq_id: residue.id().0,
        label_comp_id: residue.name().unwrap_or_default().to_string(),
        label_comp_type: component.component_type().code(),
        pdbx_pdb_ins_code: insertion_code(residue),
    }
}

fn ring_ring_label(kind: RingRingInteraction) -> &'static str {
    match kind {
        RingRingInteraction::Ff => "FF",
        RingRingInteraction::Of => "OF",
        RingRingInteraction::Ee => "EE",
        RingRingInteraction::Ft => "FT",
        RingRingInteraction::Ot => "OT",
        RingRingInteraction::Et => "ET",
        RingRingInteraction::Fe => "FE",
        RingRingInteraction::Oe => "OE",
        RingRingInteraction::Ef => "EF",
    }
}

fn ring_atom_label(kind: RingAtomInteraction) -> &'static str {
    match kind {
        RingAtomInteraction::CarbonPi => "CARBONPI",
        RingAtomInteraction::CationPi => "CATIONPI",
        RingAtomInteraction::DonorPi => "DONORPI",
        RingAtomInteraction::HalogenPi => "HALOGENPI",
        RingAtomInteraction::MetSulphurPi => "METSULPHURPI",
    }
}

fn amide_member_ids(amide: &AmideGroup) -> [String; 4] {
    [
        amide.nitrogen_id.clone(),
        amide.carbon_id.clone(),
        amide.oxygen_id.clone(),
        amide.other_carbon_id.clone(),
    ]
}

// Every real, distinct residue instance in `pdb`, each carrying its chain
// (needed for `auth_asym_id`, which a bare `&Residue` doesn't have) --
// deduped by residue identity (pointer equality into the same live `PDB`),
// same technique `contacts::find_contacts` already uses for intra-residue
// filtering.
fn residue_instances(pdb: &PDB) -> Vec<AtomConformerResidueChainModel<'_>> {
    let mut seen: Vec<*const Residue> = Vec::new();
    let mut result = Vec::new();
    for hierarchy in pdb.atoms_with_hierarchy() {
        let ptr = hierarchy.residue() as *const Residue;
        if seen.contains(&ptr) {
            continue;
        }
        seen.push(ptr);
        result.push(hierarchy);
    }
    result
}

// One perceived ring, in one real residue instance, with its real
// geometry already resolved.
struct RingInstance<'a> {
    hierarchy: AtomConformerResidueChainModel<'a>,
    component: &'a CcdComponent,
    ring: RingAtoms,
    geometry: RingGeometry,
}

// Every real ring instance across the whole structure. A residue whose
// comp_id has no matching `CcdComponent` simply contributes no rings --
// unlike `export_atom_atom_contacts`, this doesn't fail loudly for that,
// since an unknown residue can't be typed/perceived at all regardless (see
// `export_atom_plane_contacts`'s own doc for why atom-plane's per-atom
// lookups follow the same "skip, don't fail" posture).
fn collect_ring_instances<'a>(
    pdb: &'a PDB,
    components: &'a HashMap<String, CcdComponent>,
) -> Vec<RingInstance<'a>> {
    let mut instances = Vec::new();
    for hierarchy in residue_instances(pdb) {
        let comp_id = hierarchy.residue().name().unwrap_or_default();
        let Some(component) = components.get(comp_id) else {
            continue;
        };
        for ring in rings::perceive_rings(component) {
            if let Some(geometry) = rings::ring_geometry(&ring, hierarchy.residue()) {
                instances.push(RingInstance {
                    hierarchy: hierarchy.clone(),
                    component,
                    ring,
                    geometry,
                });
            }
        }
    }
    instances
}

struct AmideInstance<'a> {
    hierarchy: AtomConformerResidueChainModel<'a>,
    component: &'a CcdComponent,
    amide: AmideGroup,
    geometry: AmideGeometry,
}

fn collect_amide_instances<'a>(
    pdb: &'a PDB,
    components: &'a HashMap<String, CcdComponent>,
) -> Vec<AmideInstance<'a>> {
    let mut instances = Vec::new();
    for hierarchy in residue_instances(pdb) {
        let comp_id = hierarchy.residue().name().unwrap_or_default();
        let Some(component) = components.get(comp_id) else {
            continue;
        };
        for amide in rings::perceive_amide_groups(component) {
            if let Some(geometry) = rings::amide_geometry(&amide, hierarchy.residue()) {
                instances.push(AmideInstance {
                    hierarchy: hierarchy.clone(),
                    component,
                    amide,
                    geometry,
                });
            }
        }
    }
    instances
}

// Real pdbe-arpeggio's plane-plane (ring-ring, `interactions.py`'s
// `__calculate_plane_plane_contacts`): every pair of distinct real ring
// instances, both sides in `selection_plus` (real arpeggio's own
// `ring_key not in self.selection_plus_ring_ids or ring_key2 not in ...`
// filter -- a ring counts as "in" if any real atom of its owning residue
// is), within `rings::CENTROID_DISTANCE_MAX`, excluding an intra-residue
// `EE` pair specifically (real arpeggio's own "don't count intra-residue
// edge-to-edge to avoid intra-heterocycle interactions" rule -- everything
// else, including other intra-residue pairs, is kept, matching the
// oracle).
pub fn export_plane_plane_contacts(
    pdb: &PDB,
    components: &HashMap<String, CcdComponent>,
    selection: &crate::selection::SelectionContext,
) -> Vec<PlanePlaneContactJson> {
    let ring_instances: Vec<_> = collect_ring_instances(pdb, components)
        .into_iter()
        .filter(|r| {
            selection
                .selection_plus
                .contains_any_atom_of(r.hierarchy.residue())
        })
        .collect();
    // M8: the O(n^2) ring-pair scan is the expensive part here (real
    // structures can have hundreds of real rings once nucleotide bases and
    // ligands are counted) -- each pair's classification is fully
    // independent of every other, so the outer index runs across a rayon
    // thread pool, same drop-in-parallel shape as `find_contacts`.
    (0..ring_instances.len())
        .into_par_iter()
        .flat_map_iter(|i| {
            let ring_instances = &ring_instances;
            let selection = &selection;
            ((i + 1)..ring_instances.len()).filter_map(move |j| {
                let a = &ring_instances[i];
                let b = &ring_instances[j];
                let (distance, kind) = rings::classify_ring_ring(&a.geometry, &b.geometry)?;
                let intra_residue = std::ptr::eq(a.hierarchy.residue(), b.hierarchy.residue());
                if intra_residue && kind == RingRingInteraction::Ee {
                    return None;
                }

                let interacting_entities = interacting_entities_group(
                    selection
                        .selection
                        .contains_any_atom_of(a.hierarchy.residue()),
                    selection
                        .selection
                        .contains_any_atom_of(b.hierarchy.residue()),
                    true,
                    true,
                );

                Some(PlanePlaneContactJson {
                    bgn: group_identity(&a.hierarchy, a.component, &a.ring.atom_ids),
                    end: group_identity(&b.hierarchy, b.component, &b.ring.atom_ids),
                    entry_type: "plane-plane",
                    distance: round_2(distance),
                    contact: vec![ring_ring_label(kind)],
                    interacting_entities,
                })
            })
        })
        .collect()
}

// Real pdbe-arpeggio's atom-plane (ring-atom, `__calculate_atom_plane_contacts`):
// every real (non-hydrogen) atom in `selection_plus` against every real
// ring instance in `selection_plus`. Unlike `export_atom_atom_contacts`,
// an atom whose own residue isn't known is silently skipped rather than
// failing the whole export: real arpeggio's own SIFt-typing step can't
// classify it either way (there's no chemistry to check the ring-face
// angle/typing rules against), so skipping it loses nothing an error
// would have preserved -- it simply never could have produced an
// interaction.
pub fn export_atom_plane_contacts(
    pdb: &PDB,
    components: &HashMap<String, CcdComponent>,
    selection: &crate::selection::SelectionContext,
) -> Vec<AtomPlaneContactJson> {
    let ring_instances: Vec<_> = collect_ring_instances(pdb, components)
        .into_iter()
        .filter(|r| {
            selection
                .selection_plus
                .contains_any_atom_of(r.hierarchy.residue())
        })
        .collect();
    // M8: parallelized over real atoms rather than rings -- there are
    // typically far more real atoms than real rings in a structure, so
    // this axis balances better across a thread pool. `join_atom`/
    // `type_atom` don't depend on which ring is being checked, so they're
    // computed once per atom up front (also fixes a real inefficiency the
    // original nested-loop version had: calling them once per (ring,
    // atom) pair instead of once per atom).
    let atoms: Vec<_> = pdb.atoms_with_hierarchy().collect();

    atoms
        .into_par_iter()
        .filter(|h| h.atom().element().map(|e| e.symbol()) != Some("H"))
        .filter(|h| selection.selection_plus.contains(h.atom()))
        .filter_map(|hierarchy| {
            let joined = crate::join::join_atom(hierarchy.clone(), components);
            let (Some(ccd_atom), Some(atom_component)) = (joined.ccd_atom, joined.component) else {
                return None;
            };
            Some((hierarchy, ccd_atom, atom_component))
        })
        .flat_map_iter(|(hierarchy, ccd_atom, atom_component)| {
            let ring_instances = &ring_instances;
            let selection = &selection;
            let bits = crate::typing::type_atom(ccd_atom, atom_component);
            ring_instances.iter().filter_map(move |ring| {
                let interactions = rings::classify_ring_atom(
                    &ring.geometry,
                    hierarchy.atom().pos(),
                    ccd_atom.element(),
                    hierarchy.residue().name().unwrap_or_default(),
                    bits,
                );
                if interactions.is_empty() {
                    return None;
                }
                let mut labels: Vec<&'static str> =
                    interactions.iter().map(|i| ring_atom_label(*i)).collect();
                labels.sort_unstable();

                let distance = euclidean_distance(hierarchy.atom().pos(), ring.geometry.center);

                let interacting_entities = interacting_entities_group(
                    selection.selection.contains(hierarchy.atom()),
                    selection
                        .selection
                        .contains_any_atom_of(ring.hierarchy.residue()),
                    true,
                    true,
                );

                Some(AtomPlaneContactJson {
                    bgn: atom_identity(&hierarchy, atom_component),
                    end: group_identity(&ring.hierarchy, ring.component, &ring.ring.atom_ids),
                    entry_type: "atom-plane",
                    distance: round_2(distance),
                    contact: labels,
                    interacting_entities,
                })
            })
        })
        .collect()
}

// Real pdbe-arpeggio's group-group (amide-amide, `__calculate_group_group_contacts`):
// every pair of distinct real amide instances, both in `selection_plus`,
// that pass the shared face-on geometric test.
pub fn export_group_group_contacts(
    pdb: &PDB,
    components: &HashMap<String, CcdComponent>,
    selection: &crate::selection::SelectionContext,
) -> Vec<PlanePlaneContactJson> {
    let amide_instances: Vec<_> = collect_amide_instances(pdb, components)
        .into_iter()
        .filter(|a| {
            selection
                .selection_plus
                .contains_any_atom_of(a.hierarchy.residue())
        })
        .collect();
    // M8: same drop-in parallel shape as `export_plane_plane_contacts` --
    // amide instances are typically few, so this particular O(n^2) scan
    // is rarely a real hotspot, but the pattern stays uniform across every
    // pairwise export function rather than parallelizing some and not
    // others.
    (0..amide_instances.len())
        .into_par_iter()
        .flat_map_iter(|i| {
            let amide_instances = &amide_instances;
            let selection = &selection;
            ((i + 1)..amide_instances.len()).filter_map(move |j| {
                let a = &amide_instances[i];
                let b = &amide_instances[j];
                let distance = rings::classify_amide_amide(&a.geometry, &b.geometry)?;
                let interacting_entities = interacting_entities_group(
                    selection
                        .selection
                        .contains_any_atom_of(a.hierarchy.residue()),
                    selection
                        .selection
                        .contains_any_atom_of(b.hierarchy.residue()),
                    true,
                    true,
                );
                Some(PlanePlaneContactJson {
                    bgn: group_identity(&a.hierarchy, a.component, &amide_member_ids(&a.amide)),
                    end: group_identity(&b.hierarchy, b.component, &amide_member_ids(&b.amide)),
                    entry_type: "group-group",
                    distance: round_2(distance),
                    contact: vec!["AMIDEAMIDE"],
                    interacting_entities,
                })
            })
        })
        .collect()
}

// Real pdbe-arpeggio's group-plane (amide-ring, `__calculate_group_plane_contacts`):
// every real amide instance in `selection_plus` against every real ring
// instance in `selection_plus` that passes the shared face-on geometric
// test.
pub fn export_group_plane_contacts(
    pdb: &PDB,
    components: &HashMap<String, CcdComponent>,
    selection: &crate::selection::SelectionContext,
) -> Vec<PlanePlaneContactJson> {
    let amide_instances: Vec<_> = collect_amide_instances(pdb, components)
        .into_iter()
        .filter(|a| {
            selection
                .selection_plus
                .contains_any_atom_of(a.hierarchy.residue())
        })
        .collect();
    let ring_instances: Vec<_> = collect_ring_instances(pdb, components)
        .into_iter()
        .filter(|r| {
            selection
                .selection_plus
                .contains_any_atom_of(r.hierarchy.residue())
        })
        .collect();
    // M8: parallelized over amides (usually the smaller of the two axes)
    // with rings inner.
    amide_instances
        .par_iter()
        .flat_map_iter(|amide| {
            let ring_instances = &ring_instances;
            let selection = &selection;
            ring_instances.iter().filter_map(move |ring| {
                let distance = rings::classify_amide_ring(&amide.geometry, &ring.geometry)?;
                let interacting_entities = interacting_entities_group(
                    selection
                        .selection
                        .contains_any_atom_of(amide.hierarchy.residue()),
                    selection
                        .selection
                        .contains_any_atom_of(ring.hierarchy.residue()),
                    true,
                    true,
                );
                Some(PlanePlaneContactJson {
                    bgn: group_identity(
                        &amide.hierarchy,
                        amide.component,
                        &amide_member_ids(&amide.amide),
                    ),
                    end: group_identity(&ring.hierarchy, ring.component, &ring.ring.atom_ids),
                    entry_type: "group-plane",
                    distance: round_2(distance),
                    contact: vec!["AMIDERING"],
                    interacting_entities,
                })
            })
        })
        .collect()
}

// Every contact of every type this module can export, combined into one
// list the way real `Arpeggio.get_contacts` returns them all together.
// `Err` only from the atom-atom half (see `export_atom_atom_contacts`);
// the plane/group halves never fail loudly (see their own docs).
#[derive(Debug, Clone, Serialize)]
#[serde(untagged)]
pub enum ContactJson {
    AtomAtom(AtomAtomContactJson),
    Group(PlanePlaneContactJson),
    AtomPlane(AtomPlaneContactJson),
}

pub fn export_all_contacts(
    pdb: &PDB,
    components: &HashMap<String, CcdComponent>,
    selection: &crate::selection::SelectionContext,
) -> Result<Vec<ContactJson>, String> {
    let mut entries: Vec<ContactJson> = export_atom_atom_contacts(pdb, components, selection)?
        .into_iter()
        .map(ContactJson::AtomAtom)
        .collect();
    entries.extend(
        export_plane_plane_contacts(pdb, components, selection)
            .into_iter()
            .map(ContactJson::Group),
    );
    entries.extend(
        export_atom_plane_contacts(pdb, components, selection)
            .into_iter()
            .map(ContactJson::AtomPlane),
    );
    entries.extend(
        export_group_group_contacts(pdb, components, selection)
            .into_iter()
            .map(ContactJson::Group),
    );
    entries.extend(
        export_group_plane_contacts(pdb, components, selection)
            .into_iter()
            .map(ContactJson::Group),
    );
    Ok(entries)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rspeggio_ccd::common::common_components;

    #[test]
    fn a_real_1ca2_contact_matches_the_golden_fixtures_bgn_end_distance_and_contact_list() {
        // VAL143 CG1 <-> VAL121 CG2, taken directly from
        // tests/fixtures/golden/1CA2.json: {"bgn": {CG1, VAL 143}, "end":
        // {CG2, VAL 121}, "distance": 4.38, "contact": ["proximal",
        // "hydrophobic"]}. Deliberately picked a HYDROPHOBIC (distance +
        // typing only) pair rather than an HBOND/WEAK_HBOND one: those
        // depend on this project's own analytically-*placed* donor
        // hydrogens (`hydrogenate.rs`), which real arpeggio instead gets
        // from OpenBabel's own placement algorithm -- confirmed directly
        // (via an ad-hoc angle check before writing this test) that a
        // nearby real GLU117 O / PHE95 CA pair is exactly this kind of
        // near-threshold case: this project's placed HA gives a WEAK_HBOND
        // angle of 133.3 degrees (passes our 130 degree cutoff), while the
        // golden fixture's real OpenBabel-derived placement apparently
        // didn't. A principled, expected divergence (see HANDOFF.md), not
        // a bug -- but the wrong case to assert byte-parity on.
        // `interacting_entities` there is "INTRA_NON_SELECTION" (that
        // fixture's real, unreproduced selection) -- not checked here, see
        // module doc.
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1CA2.cif").expect("1CA2 should load");
        let mut components = common_components();
        let zn =
            rspeggio_ccd::parser::load_ccd_component("../rspeggio-ccd/tests/fixtures/ZN_ideal.cif")
                .expect("ZN fixture should parse");
        components.insert("ZN".to_string(), zn);
        let selection = crate::selection::SelectionContext::whole_structure(&pdb);

        let entries = export_atom_atom_contacts(&pdb, &components, &selection)
            .expect("1CA2 is amino acids + water + zinc, all now known");

        let entry = entries
            .iter()
            .find(|e| {
                (e.bgn.label_comp_id == "VAL"
                    && e.bgn.auth_seq_id == 143
                    && e.bgn.auth_atom_id == "CG1"
                    && e.end.label_comp_id == "VAL"
                    && e.end.auth_seq_id == 121
                    && e.end.auth_atom_id == "CG2")
                    || (e.end.label_comp_id == "VAL"
                        && e.end.auth_seq_id == 143
                        && e.end.auth_atom_id == "CG1"
                        && e.bgn.label_comp_id == "VAL"
                        && e.bgn.auth_seq_id == 121
                        && e.bgn.auth_atom_id == "CG2")
            })
            .expect("VAL143 CG1 <-> VAL121 CG2 should be a real contact in 1CA2");

        assert_eq!(entry.distance, 4.38);
        assert_eq!(entry.contact, vec!["proximal", "hydrophobic"]);
    }

    #[test]
    fn a_real_1mbo_atom_identity_matches_the_golden_fixture() {
        // VAL68 CG1, taken from the atom-plane entry in
        // tests/fixtures/golden/1MBO.json used as the worked example in
        // `rings.rs`'s own design discussion -- confirms `atom_identity`'s
        // field construction (chain, seq id, comp id/type, ins code)
        // independent of the atom-atom-specific contact list.
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1MBO.cif").expect("1MBO should load");
        let mut components = common_components();
        for (comp_id, path) in [
            ("HEM", "tests/fixtures/ccd/HEM.cif"),
            ("OXY", "tests/fixtures/ccd/OXY.cif"),
        ] {
            let component = rspeggio_ccd::parser::load_ccd_component(path)
                .unwrap_or_else(|| panic!("{comp_id} fixture should parse"));
            components.insert(comp_id.to_string(), component);
        }
        let selection = crate::selection::SelectionContext::whole_structure(&pdb);

        let entries = export_atom_atom_contacts(&pdb, &components, &selection)
            .expect("1MBO is amino acids + water + heme + bound oxygen, all now known");

        let val68_cg1 = entries
            .iter()
            .flat_map(|e| [&e.bgn, &e.end])
            .find(|a| a.label_comp_id == "VAL" && a.auth_seq_id == 68 && a.auth_atom_id == "CG1")
            .expect("VAL68 CG1 should appear in at least one real contact");

        assert_eq!(val68_cg1.auth_asym_id, "A");
        assert_eq!(val68_cg1.label_comp_type, "P");
        assert_eq!(val68_cg1.pdbx_pdb_ins_code, " ");
    }

    #[test]
    fn a_real_water_water_contact_is_classified_correctly_in_whole_structure_mode() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");
        let components = common_components();
        let selection = crate::selection::SelectionContext::whole_structure(&pdb);

        let entries =
            export_atom_atom_contacts(&pdb, &components, &selection).expect("1UBQ should export");

        let found = entries.iter().any(|e| {
            e.bgn.label_comp_id == "HOH"
                && e.end.label_comp_id == "HOH"
                && e.interacting_entities == "WATER_WATER"
        });
        assert!(
            found,
            "expected at least one real water-water contact in 1UBQ"
        );
    }

    #[test]
    fn a_real_protein_water_contact_is_selection_water_in_whole_structure_mode() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");
        let components = common_components();
        let selection = crate::selection::SelectionContext::whole_structure(&pdb);

        let entries =
            export_atom_atom_contacts(&pdb, &components, &selection).expect("1UBQ should export");

        let found = entries.iter().any(|e| {
            (e.bgn.label_comp_id == "HOH") != (e.end.label_comp_id == "HOH")
                && e.interacting_entities == "SELECTION_WATER"
        });
        assert!(
            found,
            "expected at least one real protein-water contact in 1UBQ"
        );
    }

    #[test]
    fn a_real_protein_protein_contact_is_intra_selection_in_whole_structure_mode() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");
        let components = common_components();
        let selection = crate::selection::SelectionContext::whole_structure(&pdb);

        let entries =
            export_atom_atom_contacts(&pdb, &components, &selection).expect("1UBQ should export");

        let found = entries.iter().any(|e| {
            e.bgn.label_comp_id != "HOH"
                && e.end.label_comp_id != "HOH"
                && e.interacting_entities == "INTRA_SELECTION"
        });
        assert!(
            found,
            "expected at least one real protein-protein contact in 1UBQ"
        );
    }

    #[test]
    fn every_entry_has_exactly_one_distance_category_label() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");
        let components = common_components();
        let selection = crate::selection::SelectionContext::whole_structure(&pdb);

        let entries =
            export_atom_atom_contacts(&pdb, &components, &selection).expect("1UBQ should export");
        assert!(!entries.is_empty());

        const DISTANCE_LABELS: [&str; 5] = ["clash", "covalent", "vdw_clash", "vdw", "proximal"];
        for entry in &entries {
            let category_count = entry
                .contact
                .iter()
                .filter(|label| DISTANCE_LABELS.contains(label))
                .count();
            assert_eq!(
                category_count, 1,
                "entry {entry:?} should have exactly one distance-category label"
            );
        }
    }

    #[test]
    fn an_unknown_component_fails_the_whole_export_loudly() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");
        // Deliberately empty -- simulates every residue being unknown, the
        // same way `join.rs`'s own equivalent test does.
        let components: HashMap<String, CcdComponent> = HashMap::new();
        let selection = crate::selection::SelectionContext::whole_structure(&pdb);

        assert!(export_atom_atom_contacts(&pdb, &components, &selection).is_err());
    }

    #[test]
    fn the_exported_shape_actually_serializes_to_the_real_json_field_names() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1UBQ.cif").expect("1UBQ should load");
        let components = common_components();
        let selection = crate::selection::SelectionContext::whole_structure(&pdb);

        let entries =
            export_atom_atom_contacts(&pdb, &components, &selection).expect("1UBQ should export");
        let json = serde_json::to_value(&entries[0]).expect("should serialize");

        for field in [
            "bgn",
            "end",
            "type",
            "distance",
            "contact",
            "interacting_entities",
        ] {
            assert!(
                json.get(field).is_some(),
                "expected real pdbe-arpeggio field {field:?} in serialized output"
            );
        }
        assert!(json["bgn"].get("pdbx_PDB_ins_code").is_some());
        assert_eq!(json["type"], "atom-atom");
    }

    #[test]
    fn a_real_ring_ring_contact_matches_rings_rs_own_worked_example() {
        // TRP5's 5-membered pyrrole ring <-> TRP16's own, EF, 5.58A -- the
        // exact real pair `rings.rs`'s ring-perception/geometry work
        // already surfaced from 1CA2 (see this crate's own commit
        // history), re-verified here through the export path instead of
        // calling `rings::classify_ring_ring` directly.
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1CA2.cif").expect("1CA2 should load");
        let components = common_components();
        let selection = crate::selection::SelectionContext::whole_structure(&pdb);

        let entries = export_plane_plane_contacts(&pdb, &components, &selection);
        let entry = entries
            .iter()
            .find(|e| {
                e.bgn.label_comp_id == "TRP"
                    && e.bgn.auth_seq_id == 5
                    && e.end.label_comp_id == "TRP"
                    && e.end.auth_seq_id == 16
                    && e.bgn.auth_atom_id == "CD1,CD2,CE2,CG,NE1"
                    && e.end.auth_atom_id == "CD1,CD2,CE2,CG,NE1"
            })
            .expect("TRP5's pyrrole ring <-> TRP16's pyrrole ring should be a real contact");

        assert_eq!(entry.entry_type, "plane-plane");
        assert_eq!(entry.distance, 5.58);
        assert_eq!(entry.contact, vec!["EF"]);
        assert_eq!(entry.interacting_entities, "INTRA_SELECTION");
    }

    #[test]
    fn a_real_amide_amide_contact_matches_rings_rs_own_worked_example() {
        // ASN67 <-> GLN92, AMIDEAMIDE, 4.23A -- same real pair `rings.rs`'s
        // own amide-amide test uses.
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1CA2.cif").expect("1CA2 should load");
        let components = common_components();
        let selection = crate::selection::SelectionContext::whole_structure(&pdb);

        let entries = export_group_group_contacts(&pdb, &components, &selection);
        let entry = entries
            .iter()
            .find(|e| {
                e.bgn.label_comp_id == "ASN"
                    && e.bgn.auth_seq_id == 67
                    && e.end.label_comp_id == "GLN"
                    && e.end.auth_seq_id == 92
            })
            .expect("ASN67 <-> GLN92 should be a real amide-amide contact");

        assert_eq!(entry.entry_type, "group-group");
        assert_eq!(entry.distance, 4.23);
        assert_eq!(entry.contact, vec!["AMIDEAMIDE"]);
        assert_eq!(entry.bgn.auth_atom_id, "CB,CG,ND2,OD1");
        assert_eq!(entry.end.auth_atom_id, "CD,CG,NE2,OE1");
    }

    #[test]
    fn a_real_carbon_pi_atom_plane_contact_is_found_in_1ca2() {
        // PRO201's ring atoms sitting close to TRP5's pyrrole face --
        // confirmed via an ad-hoc probe before writing this test (83 real
        // atom-plane entries in 1CA2 alone; this is the first).
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1CA2.cif").expect("1CA2 should load");
        let components = common_components();
        let selection = crate::selection::SelectionContext::whole_structure(&pdb);

        let entries = export_atom_plane_contacts(&pdb, &components, &selection);
        let entry = entries
            .iter()
            .find(|e| {
                e.bgn.label_comp_id == "PRO"
                    && e.bgn.auth_seq_id == 201
                    && e.bgn.auth_atom_id == "CB"
                    && e.end.label_comp_id == "TRP"
                    && e.end.auth_seq_id == 5
                    && e.end.auth_atom_id == "CD1,CD2,CE2,CG,NE1"
            })
            .expect("PRO201 CB <-> TRP5's pyrrole ring should be a real atom-plane contact");

        assert_eq!(entry.entry_type, "atom-plane");
        assert_eq!(entry.distance, 3.99);
        assert_eq!(entry.contact, vec!["CARBONPI"]);
    }

    #[test]
    fn hem_only_contributes_its_two_flagged_pyrrole_rings_a_real_documented_limitation() {
        // See this module's own doc comment: HEM's real CCD entry flags
        // aromaticity inconsistently across its own four chemically
        // equivalent pyrrole rings -- rings A and C are flagged, B and D
        // aren't -- so this project's flag-based ring perception sees only
        // 2 of HEM's 4 real rings, unlike real pdbe-arpeggio's
        // OpenBabel-based re-perception (which golden fixture 1MBO.json's
        // own real CARBONPI-against-ring-B entry confirms treats all four
        // consistently). This test pins the real, precise shape of that
        // gap down rather than leaving it as an unverified claim in a
        // comment.
        let hem_component = rspeggio_ccd::parser::load_ccd_component("tests/fixtures/ccd/HEM.cif")
            .expect("HEM fixture should parse");

        let mut ring_atom_sets: Vec<Vec<String>> = rings::perceive_rings(&hem_component)
            .into_iter()
            .map(|r| {
                let mut ids = r.atom_ids;
                ids.sort();
                ids
            })
            .collect();
        ring_atom_sets.sort();

        assert_eq!(
            ring_atom_sets,
            vec![
                vec!["C1A", "C2A", "C3A", "C4A", "NA"],
                vec!["C1C", "C2C", "C3C", "C4C", "NC"],
            ],
            "expected exactly HEM's flagged A and C pyrrole rings, not its unflagged B/D ones"
        );
    }

    #[test]
    fn export_all_contacts_combines_every_real_contact_type() {
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1CA2.cif").expect("1CA2 should load");
        let mut components = common_components();
        let zn =
            rspeggio_ccd::parser::load_ccd_component("../rspeggio-ccd/tests/fixtures/ZN_ideal.cif")
                .expect("ZN fixture should parse");
        components.insert("ZN".to_string(), zn);
        let selection = crate::selection::SelectionContext::whole_structure(&pdb);

        let entries =
            export_all_contacts(&pdb, &components, &selection).expect("1CA2 should export");

        let mut seen_types: Vec<&str> = entries
            .iter()
            .map(|e| match e {
                ContactJson::AtomAtom(_) => "atom-atom",
                ContactJson::Group(g) => g.entry_type,
                ContactJson::AtomPlane(_) => "atom-plane",
            })
            .collect();
        seen_types.sort_unstable();
        seen_types.dedup();

        // group-plane happens not to occur anywhere in 1CA2 (confirmed via
        // the same ad-hoc probe as the other tests above), so it's
        // deliberately not asserted for here.
        for expected in ["atom-atom", "plane-plane", "atom-plane", "group-group"] {
            assert!(
                seen_types.contains(&expected),
                "expected at least one real {expected} contact in the combined export, got types {seen_types:?}"
            );
        }
    }

    // The tests below use a real `RESNAME:` selection instead of
    // whole-structure mode -- and so are checkable against
    // `tests/fixtures/golden/*.json`'s own real `interacting_entities`
    // values directly, not just schema/bgn/end/distance/contact (as the
    // whole-structure-mode tests above are limited to). These 4 selections
    // (ZN/HEM/FMN/FES) aren't guessed: each was worked out by hand-checking
    // which single real comp_id, when treated as "the selection", makes
    // every golden `interacting_entities` value in that fixture consistent
    // (e.g. 1CA2: every ZN-involving entry is INTER or SELECTION_WATER,
    // every non-ZN entry is INTRA_NON_SELECTION or NON_SELECTION_WATER,
    // never INTRA_SELECTION since there's only one real zinc) -- see this
    // module's own doc comment for the reasoning.

    #[test]
    fn a_real_zn_selection_reproduces_1ca2s_golden_interacting_entities_exactly() {
        // GLU117 O <-> PHE95 CA is real pdbe-arpeggio's own
        // INTRA_NON_SELECTION example (neither atom near the zinc site);
        // ZN <-> HIS96 ND1 is a real INTER example (the zinc-coordinating
        // histidine). Both confirmed directly against
        // tests/fixtures/golden/1CA2.json's own real values -- true
        // byte-level parity on `interacting_entities`, not just a
        // plausible-looking label.
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1CA2.cif").expect("1CA2 should load");
        let mut components = common_components();
        let zn =
            rspeggio_ccd::parser::load_ccd_component("../rspeggio-ccd/tests/fixtures/ZN_ideal.cif")
                .expect("ZN fixture should parse");
        components.insert("ZN".to_string(), zn);

        let selection =
            crate::selection::SelectionContext::from_specs(&pdb, &["RESNAME:ZN".to_string()])
                .expect("ZN should be found in 1CA2");
        let entries = export_atom_atom_contacts(&pdb, &components, &selection)
            .expect("1CA2 is amino acids + water + zinc, all now known");

        let glu_phe = entries
            .iter()
            .find(|e| {
                (e.bgn.label_comp_id == "GLU"
                    && e.bgn.auth_seq_id == 117
                    && e.bgn.auth_atom_id == "O"
                    && e.end.label_comp_id == "PHE"
                    && e.end.auth_seq_id == 95
                    && e.end.auth_atom_id == "CA")
                    || (e.end.label_comp_id == "GLU"
                        && e.end.auth_seq_id == 117
                        && e.end.auth_atom_id == "O"
                        && e.bgn.label_comp_id == "PHE"
                        && e.bgn.auth_seq_id == 95
                        && e.bgn.auth_atom_id == "CA")
            })
            .expect("GLU117 O <-> PHE95 CA should be a real contact in 1CA2");
        assert_eq!(glu_phe.interacting_entities, "INTRA_NON_SELECTION");

        let zn_his = entries
            .iter()
            .find(|e| {
                (e.bgn.label_comp_id == "ZN"
                    && e.end.label_comp_id == "HIS"
                    && e.end.auth_seq_id == 96
                    && e.end.auth_atom_id == "ND1")
                    || (e.end.label_comp_id == "ZN"
                        && e.bgn.label_comp_id == "HIS"
                        && e.bgn.auth_seq_id == 96
                        && e.bgn.auth_atom_id == "ND1")
            })
            .expect("ZN <-> HIS96 ND1 should be a real contact in 1CA2");
        assert_eq!(zn_his.interacting_entities, "INTER");
    }

    #[test]
    fn a_real_hem_selection_reproduces_1mbos_golden_interacting_entities_exactly() {
        // HEM's own real bound OXY is INTER (only HEM is selected, not
        // OXY) in the golden fixture -- confirmed directly against
        // tests/fixtures/golden/1MBO.json (every real HEM-OXY pair there
        // is INTER, never INTRA_SELECTION, since only one of the two is
        // ever selected).
        let (pdb, _errors) =
            pdbtbx::open("../../tests/fixtures/structures/1MBO.cif").expect("1MBO should load");
        let mut components = common_components();
        for (comp_id, path) in [
            ("HEM", "tests/fixtures/ccd/HEM.cif"),
            ("OXY", "tests/fixtures/ccd/OXY.cif"),
        ] {
            let component = rspeggio_ccd::parser::load_ccd_component(path)
                .unwrap_or_else(|| panic!("{comp_id} fixture should parse"));
            components.insert(comp_id.to_string(), component);
        }

        let selection =
            crate::selection::SelectionContext::from_specs(&pdb, &["RESNAME:HEM".to_string()])
                .expect("HEM should be found in 1MBO");
        let entries = export_atom_atom_contacts(&pdb, &components, &selection)
            .expect("1MBO is amino acids + water + heme + bound oxygen, all now known");

        let found_inter = entries.iter().any(|e| {
            let ids = [e.bgn.label_comp_id.as_str(), e.end.label_comp_id.as_str()];
            ids.contains(&"HEM") && ids.contains(&"OXY") && e.interacting_entities == "INTER"
        });
        assert!(
            found_inter,
            "expected a real HEM-OXY contact classified INTER (only HEM is selected)"
        );
    }
}
