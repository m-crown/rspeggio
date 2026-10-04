// crates/rspeggio-py/src/lib.rs
//
// M9: PyO3 bindings -- exposes this project's contact-detection pipeline
// (M5-M8) to Python as a native extension module.
//
// Scope decision: a single high-level function, `get_contacts`, returning
// the same real pdbe-arpeggio-shaped JSON `export.rs` already produces
// (`export::export_all_contacts`), as a JSON *string* rather than parsed
// Python objects. This avoids needing a Rust<->Python value-bridging
// dependency (e.g. `pythonize`) for what's currently a single call --
// `json.loads()` on the Python side is one line, and the real schema
// (field names, nesting) is already exactly what real pdbe-arpeggio's own
// JSON output looks like, so there's nothing Python-specific to add on
// top of it. A richer typed Python API (dataclasses per contact type,
// etc.) is real, deferred scope, not attempted here.
//
// Structure/CCD loading happens on this side (not in `rspeggio-core`,
// which stays a pure library with no I/O policy of its own) -- same
// `common_components()` bundled set every Rust-side test already uses,
// plus caller-supplied extra CCD files for any hetero groups/ligands not
// in it (decision 03: unknown components fail loudly, surfaced here as a
// real Python exception rather than a silent gap). `geometric_fallback`
// opts out of that for anything still unknown after loading, via
// `rspeggio_core::perception` (M10) -- off by default, so the fail-loudly
// default is unchanged.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rspeggio_ccd::common::common_components;
use rspeggio_ccd::component::CcdComponent;
use rspeggio_core::export::export_all_contacts;
use rspeggio_core::selection::SelectionContext;
use std::collections::HashMap;

// Loads the bundled common CCD set, then merges in any caller-supplied
// extra CCD `.cif` files (e.g. a bound ligand's real CCD entry), keyed by
// each file's own real comp_id (`_chem_comp.id`) -- not its filename.
// Real CCD downloads don't reliably name the file after the comp_id (RCSB's
// own "ideal coordinates" convention names Zn2+'s entry `ZN_ideal.cif`,
// whose real id is `"ZN"`), confirmed directly: an earlier filename-stem
// version of this function broke on exactly that file during this
// binding's own end-to-end testing.
fn load_components(extra_ccd_paths: &[String]) -> PyResult<HashMap<String, CcdComponent>> {
    let mut components = common_components();
    for path in extra_ccd_paths {
        let (comp_id, component) = rspeggio_ccd::parser::load_ccd_component_with_id(path)
            .ok_or_else(|| PyValueError::new_err(format!("failed to parse CCD file: {path}")))?;
        components.insert(comp_id, component);
    }
    Ok(components)
}

/// Computes every real contact in a structure, in real pdbe-arpeggio's own
/// JSON schema (atom-atom, plane-plane, atom-plane, group-group,
/// group-plane), returned as a JSON string.
///
/// Args:
///     structure_path: path to a PDB or mmCIF structure file.
///     extra_ccd_paths: paths to additional CCD component ``.cif`` files
///         for any hetero groups/ligands not already in the bundled
///         common set (the 20 standard amino acids, water, standard
///         nucleotides, and a handful of common ions). Each file is keyed
///         by its own real ``_chem_comp.id``, not its filename (RCSB's
///         ``ZN_ideal.cif`` is ``"ZN"``).
///     selection: real pdbe-arpeggio selection strings, e.g.
///         ``["RESNAME:ZN"]`` or ``["/A/45/"]``. Omit (or pass ``None``)
///         for whole-structure mode, where every residue is treated as
///         selected.
///     ph: solution pH used for the one pH-dependent typing rule this
///         project has (free/side-chain carboxyl-group protonation --
///         see ``rspeggio_core::typing::type_atom``'s own doc). Defaults
///         to physiological pH (7.4), real pdbe-arpeggio's own default.
///     geometric_fallback: if true, any residue whose comp_id still has no
///         CCD component (bundled or supplied) gets its chemistry perceived
///         from its own 3D coordinates instead of raising. Approximate --
///         no aromaticity, no charge model; see the README's
///         "Geometric-perception fallback" section for every limitation.
///         Defaults to false.
///
/// Returns:
///     A JSON string: a list of contact objects. Parse it with
///     ``json.loads`` on the Python side.
///
/// Raises:
///     ValueError: the structure file couldn't be opened, a residue's
///         comp_id has no matching CCD component (real or supplied) and
///         ``geometric_fallback`` is false, the
///         selection matched no real atoms, or a selection string is
///         malformed.
#[pyfunction]
#[pyo3(signature = (structure_path, extra_ccd_paths=Vec::new(), selection=None, ph=rspeggio_core::typing::PHYSIOLOGICAL_PH, geometric_fallback=false))]
fn get_contacts(
    structure_path: String,
    extra_ccd_paths: Vec<String>,
    selection: Option<Vec<String>>,
    ph: f64,
    geometric_fallback: bool,
) -> PyResult<String> {
    let (pdb, _errors) = pdbtbx::open(&structure_path).map_err(|errors| {
        PyValueError::new_err(format!("failed to open {structure_path:?}: {errors:?}"))
    })?;

    let mut components = load_components(&extra_ccd_paths)?;
    if geometric_fallback {
        rspeggio_core::perception::augment_with_geometric_fallback(&pdb, &mut components);
    }

    let selection_context = match selection {
        Some(specs) => SelectionContext::from_specs(&pdb, &specs).map_err(PyValueError::new_err)?,
        None => SelectionContext::whole_structure(&pdb),
    };

    let contacts = export_all_contacts(&pdb, &components, &selection_context, ph)
        .map_err(PyValueError::new_err)?;

    serde_json::to_string(&contacts)
        .map_err(|e| PyValueError::new_err(format!("failed to serialize contacts: {e}")))
}

#[pymodule]
fn rspeggio(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(get_contacts, m)?)?;
    Ok(())
}
