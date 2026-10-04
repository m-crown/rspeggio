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

`get_contacts` returns a JSON string — a list of contact objects in real pdbe-arpeggio's own schema (`atom-atom`, `plane-plane`, `atom-plane`, `group-group`, `group-plane` entry types). Errors (an unreadable structure file, a residue whose comp_id has no matching CCD component and `geometric_fallback` is off, a selection that matches no real atoms, a malformed selection string) raise a Python `ValueError`.

For a ligand you have no CCD file for, pass `geometric_fallback=True` — see [Geometric-perception fallback](#geometric-perception-fallback) for what that does and how far to trust it.

### Local Rust builds of `rspeggio-py`

PyO3's `extension-module` feature (needed to build a real, importable Python extension) deliberately does *not* link against libpython at build time — a real Python extension gets those symbols from the interpreter that loads it, not at link time. That's correct for the wheel `maturin` produces, but it means a bare `cargo build -p rspeggio-py` (or `cargo test`, or `just ci`, which touches the whole workspace) would otherwise fail with "undefined symbols" for every Python C-API call.

`.cargo/config.toml` in the repo root works around this with PyO3's own documented fix — telling the linker to defer those symbols to load time (the same thing `maturin` does under the hood) — so ordinary `cargo build`/`cargo test`/`just ci` work across the whole workspace without needing maturin installed at all. You only need maturin (and a venv) to actually produce something Python can `import`.

## Geometric-perception fallback

rspeggio gets each residue's chemistry (bonds, bond orders, aromaticity, hydrogens) from its [CCD](https://www.wwpdb.org/data/ccd) entry, looked up by comp_id. By default, a residue whose comp_id has no CCD entry — not in the bundled set, and no file supplied via `extra_ccd_paths` — **fails the whole run loudly**. That's deliberate: silently guessing chemistry is worse than an error you can fix by fetching the real CCD file (`https://files.rcsb.org/ligands/download/<ID>.cif`), which is always the better option when one exists.

For a component with no CCD entry at all (a novel or in-house ligand), you can opt in to perceiving its chemistry from its own 3D coordinates instead:

```python
contacts = json.loads(rspeggio.get_contacts("structure.cif", geometric_fallback=True))
```

```rust
let mut components = rspeggio_ccd::common::common_components();
rspeggio_core::perception::augment_with_geometric_fallback(&pdb, &mut components);
// then call the export functions as usual
```

Only comp_ids still unknown after loading are perceived; a real CCD entry is never replaced. Each perceived component is reported with `label_comp_type` `"B"` (bound molecule).

### How it works

Five passes over one residue instance's real coordinates ([`crates/rspeggio-core/src/perception.rs`](crates/rspeggio-core/src/perception.rs)):

1. **Connectivity** — two atoms in the same residue are bonded if `0.4 Å < d < r_cov₁ + r_cov₂ + 0.45 Å`, Open Babel's own `ConnectTheDots` rule (the one real pdbe-arpeggio's bond graph comes from). Bonds never cross residues.
2. **Bond order** — `Double` if `d < 0.94 × (r_cov₁ + r_cov₂)`, else `Single`. The 0.94 is calibrated against real 1UBQ bond lengths, as a ratio to the summed radii: true single bonds sit at 0.97–1.03, C=O double bonds at 0.85–0.91.
3. **Valence repair** — resonance-shortened single bonds (amide C-N 0.87–0.93, carboxylate C-O, guanidinium C-N) overlap the double-bond range, so pass 2 over-assigns doubles. While any C/N/O carries more bond order than its valence (4/3/2), its *longest* double bond is demoted to single. That recovers the one-Kekulé-form drawing a CCD entry uses: a carboxylate keeps one C=O and one C-O; an amide keeps C=O and makes C-N single.
4. **Aromaticity**, in the CCD's own convention (the one every looked-up component already follows). A 5- or 6-membered ring of C/N/O/S is aromatic if all three hold: its mean bond ratio is under 0.975, it's within 0.15 Å of flat, and no ring atom carries an exocyclic C=O or C=S. Each ring is judged on its own, as the CCD flags fused systems (guanine's 5-ring is aromatic, its 6-ring isn't). The thresholds come from real rings compared with their real CCD flags: every aromatic ring checked has a mean ratio ≤ 0.962, every saturated ring ≥ 0.991. Flatness alone can't separate them, since a real proline ring can be flatter than a real aromatic ring. The exocyclic C=O rule is what keeps uracil, thymine and cytosine non-aromatic, as in the CCD. Open Babel, and so real pdbe-arpeggio, would call them aromatic.
5. **Implicit hydrogens** — most X-ray structures contain no hydrogens, but donor typing depends on them (a CCD entry lists every H whether or not it was deposited). Each non-aromatic C/N/O/S gets `valence − bond-order sum` hydrogens. Aromatic ring atoms instead get a degree-based count: a ring C-H gets one H; a 5-ring nitrogen gets one H (pyrrole-type) unless the ring already has a substituted N, O or S; a 6-ring nitrogen gets none. The hydrogens are then placed geometrically, exactly like an undeposited CCD hydrogen. This is the role Open Babel's `AddHydrogens` plays in real pdbe-arpeggio.

### How close it gets

Measured against real CCD entries, treating known components as if unknown:

- **Ring aromaticity**: every PHE/TYR/TRP/HIS/PRO ring across the protein fixtures (165 residues) and every base in a real DNA structure (1BNA, 24 nucleotides) matches its CCD flags exactly, instance by instance. A chlorobenzene ligand (8CL) matches too.
- **Per-atom typing**, every amino-acid type in 1CA2 (127 side-chain + carbonyl-O atoms): 103 match the real CCD exactly, including every aromatic ring atom on both typing and hydrogen count. The remaining 24 are all limitations below, or cases where the CCD is the one that's off: 20 backbone O's, where the CCD's free-monomer OXT makes them look like carboxylates; LYS NZ; GLU's two equivalent carboxylate O's swapped; ARG NH2's H count. The same comparison on DNA bases matches every base atom.
- **ALA in 1UBQ**, perceived instead of looked up: the same 76 contacts, 73 labelled identically, including all hbond/polar/weak_hbond.
- **HEM in 1MBO**, a real 43-atom Fe porphyrin ligand, via the Python binding: ignoring the `aromatic` label, its output contains every contact the real CCD gives, with identical hbond, polar, ionic, metal_complex, hydrophobic, weak_hbond and weak_polar counts. Perception finds all four pyrroles aromatic, whereas HEM's CCD entry flags only two (a known curation quirk; the four are chemically equivalent). That adds 6 ring contacts and the `aromatic` label on 24 atom pairs. One of the added contacts, VAL68 CG1 → pyrrole B (CARBONPI, 3.98 Å), is the only ring contact in real pdbe-arpeggio's own 1MBO output, and it matches that entry field for field. The CCD-backed path has never produced it.

### Limitations — read before trusting the output

- **Planar conjugated rings without an exocyclic C=O can be false positives.** FMN's central N5/N10 ring is called aromatic where its CCD entry says it isn't; that's the one real false positive found.
- **The aromatic/saturated margin is narrow.** Every real ring checked falls at mean bond ratio ≤ 0.962 (aromatic) or ≥ 0.991 (saturated), with the cutoff at 0.975. A strained or low-resolution ring near the cutoff can flip.
- **Only 5- and 6-membered rings** can be aromatic, and never a ring through a metal (such as HEM's Fe chelate rings).
- **Imidazole tautomers aren't decided.** When neither ring N is substituted, both get a hydrogen (so both are donors), matching the CCD's own histidine entry.
- **No charge model.** Valences are neutral, so amines read as `-NH2`: hbond donor and acceptor, but never positively ionisable (unlike a CCD entry drawn as `-NH3+`, such as lysine's NZ). Carboxylates are unaffected, since their ionisation is decided by pH in typing.
- **Bonds to other residues aren't seen.** An atom covalently bonded to another residue (a polymer link, a covalently attached ligand) gets an implicit hydrogen in that bond's place, which can make it a spurious weak donor.
- **Amides without a third carbonyl substituent** (formamide-like) aren't repaired, so their C-N stays `Double` and the amide group isn't perceived.
- **Triple bonds** read as `Double`.
- **No over-bonding cleanup.** Open Babel follows its distance rule with a pass that removes bonds over an atom's maximum bond count or at <45° angles; that isn't ported, so a badly distorted or clashing ligand can come out over-bonded.
- **One definition per comp_id**, built from the first residue instance in the structure. Atoms missing from that instance (disorder, partial occupancy) are missing for every instance.

## Fetching new CCD/structure fixtures

Structure fixtures come straight from RCSB (`https://files.rcsb.org/download/<PDB_ID>.cif`); CCD component fixtures come from the RCSB ligand endpoint (`https://files.rcsb.org/ligands/download/<ID>.cif`). Coordinates and chemistry are always read back out of the fetched file programmatically in tests, never transcribed by hand.
