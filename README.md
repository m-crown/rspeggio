# rspeggio

A multi-threaded rust reimplementation of pdbe-arpeggio, with CLI and python bindings.

## Setup

Make sure just is installed e.g. `brew install just`

## Workspace layout

Three crates:

- **`crates/rspeggio-ccd`** — parses the Chemical Component Dictionary (CCD): atoms, bonds, and component metadata (`_chem_comp.type`/`.name`/`.id`) from real mmCIF `_chem_comp*` categories. Ships a bundled `common_components()` set (the 20 standard amino acids, water, standard nucleotides, and a handful of common ions) so the overwhelming majority of structures resolve with no network access or on-disk cache.
- **`crates/rspeggio-core`** — everything downstream of structure loading (via [`pdbtbx`](https://docs.rs/pdbtbx)) and CCD chemistry: atom typing, distance/feature contact classification, ring/amide plane geometry, the real pdbe-arpeggio JSON export, selection parsing, and unsatisfied-contact detection. This is the library; it has no I/O policy of its own beyond reading files it's given a path to.
- **`crates/rspeggio-py`** — a [PyO3](https://pyo3.rs) native extension exposing `rspeggio-core`'s pipeline to Python as a single function, `rspeggio.get_contacts(...)`. See [Python bindings](#python-bindings) below.

## Development

```bash
just check   # fmt + clippy + test — the loop to run after every change
just ci      # the exact recipe CI runs (fmt-check + clippy --all-targets -D warnings) — stricter than `check`
```

Run `just ci` (not just `just check`) before calling anything done — several real bugs have only surfaced under `-D warnings` (e.g. dead code from a function only used in tests).

Every test fixture is real data (fetched from RCSB/PDB or the real CCD, never hand-invented), and tests assert against something a self-consistency check could miss — see `tests/fixtures/` for structures and `tests/fixtures/golden/` for real pdbe-arpeggio JSON output used for parity checks.

## Python bindings

`rspeggio-py` builds a real native Python extension via [`maturin`](https://www.maturin.rs/). `cargo build`/`cargo test` work locally without maturin too (see [Local Rust builds](#local-rust-builds-of-rspeggio-py) below), but to actually use it from Python you need maturin to produce and install the wheel.

### Building and installing

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install maturin

cd crates/rspeggio-py
maturin develop        # builds the extension and installs it into the active venv, editable
# or: maturin build --release   # produces a wheel in target/wheels/ for `pip install <wheel>`
```

`maturin develop` needs an active virtualenv (`VIRTUAL_ENV` set, e.g. via `source .venv/bin/activate`) or a conda environment — it refuses to install into the system Python.

### Usage

```python
import rspeggio
import json

contacts = json.loads(rspeggio.get_contacts("structure.cif"))

# With a selection (real pdbe-arpeggio's RESNAME:<id> / /<chain>/<resnum>[<inscode>]/<atom> syntax),
# extra CCD files for hetero groups not in the bundled common set, and a non-default pH:
contacts = json.loads(rspeggio.get_contacts(
    "structure.cif",
    extra_ccd_paths=["HEM.cif"],
    selection=["RESNAME:HEM"],
    ph=7.4,
))
```

`get_contacts` returns a JSON string — a list of contact objects in real pdbe-arpeggio's own schema (`atom-atom`, `plane-plane`, `atom-plane`, `group-group`, `group-plane` entry types). Errors (an unreadable structure file, a residue whose comp_id has no matching CCD component, a selection that matches no real atoms, a malformed selection string) raise a Python `ValueError`.

### Local Rust builds of `rspeggio-py`

PyO3's `extension-module` feature (needed to build a real, importable Python extension) deliberately does *not* link against libpython at build time — a real Python extension gets those symbols from the interpreter that loads it, not at link time. That's correct for the wheel `maturin` produces, but it means a bare `cargo build -p rspeggio-py` (or `cargo test`, or `just ci`, which touches the whole workspace) would otherwise fail with "undefined symbols" for every Python C-API call.

`.cargo/config.toml` in the repo root works around this with PyO3's own documented fix — telling the linker to defer those symbols to load time (the same thing `maturin` does under the hood) — so ordinary `cargo build`/`cargo test`/`just ci` work across the whole workspace without needing maturin installed at all. You only need maturin (and a venv) to actually produce something Python can `import`.

## Fetching new CCD/structure fixtures

Structure fixtures come straight from RCSB (`https://files.rcsb.org/download/<PDB_ID>.cif`); CCD component fixtures come from the RCSB ligand endpoint (`https://files.rcsb.org/ligands/download/<ID>.cif`). Coordinates and chemistry are always read back out of the fetched file programmatically in tests, never transcribed by hand.
