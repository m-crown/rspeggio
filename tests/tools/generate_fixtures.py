"""
Regenerates tests/fixtures/golden/ by running the real `pdbe-arpeggio`
(the oracle) against each structure in tests/fixtures/structures/.

Usage (from within the oracle_env venv, with pdbe-arpeggio installed):

    python tests/tools/generate_fixtures.py

Each fixture's selection is defined in FIXTURES below - kept explicit and
version-controlled here rather than inferred, so re-running this script is
always reproducible.
"""
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
STRUCTURES_DIR = os.path.join(HERE, '..', 'fixtures', 'structures')
GOLDEN_DIR = os.path.join(HERE, '..', 'fixtures', 'golden')

# (structure filename, selection args or None for full-structure mode)
FIXTURES = [
    ('1UBQ.cif', None),               # protein only, no ligand - baseline
    ('1MBO.cif', ['RESNAME:HEM']),    # heme / porphyrin, iron
    ('1FLV.cif', ['RESNAME:FMN']),    # flavin, isoalloxazine ring
    ('4FXC.cif', ['RESNAME:FES']),    # Fe-S cluster, non-aromatic
    ('1CA2.cif', ['RESNAME:ZN']),     # single metal ion
]


def run_fixture(filename, selection):
    src = os.path.join(STRUCTURES_DIR, filename)
    if not os.path.exists(src):
        print(f'SKIP {filename}: not found in {STRUCTURES_DIR}')
        return False

    cmd = ['pdbe-arpeggio', src, '-o', GOLDEN_DIR, '-m']
    if selection:
        cmd += ['-s'] + selection

    print(f'Running: {" ".join(cmd)}')
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f'FAILED {filename}:\n{result.stdout}\n{result.stderr}')
        return False

    json_name = filename.split('.')[0] + '.json'
    out_path = os.path.join(GOLDEN_DIR, json_name)
    if not os.path.exists(out_path):
        print(f'FAILED {filename}: expected output {out_path} was not produced')
        return False

    with open(out_path) as f:
        contacts = json.load(f)
    print(f'  -> {json_name}: {len(contacts)} contact records')
    return True


def main():
    os.makedirs(GOLDEN_DIR, exist_ok=True)
    results = [run_fixture(f, sel) for f, sel in FIXTURES]
    n_ok = sum(results)
    print(f'\n{n_ok}/{len(FIXTURES)} fixtures generated successfully.')
    if n_ok != len(FIXTURES):
        sys.exit(1)


if __name__ == '__main__':
    main()
