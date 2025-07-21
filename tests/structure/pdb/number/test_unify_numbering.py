# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.


import os

from conftest import download_pdb

from simlify.structure.pdb.numbering.main import run_unify_numbering


def cache_5khb(dir_test: str) -> str:
    pdb_id = "5khb"
    path_pdb = os.path.join(dir_test, f"tmp/{pdb_id}.pdb")
    download_pdb(pdb_id, path_save=path_pdb)
    return path_pdb


def test_unify_5khb_noreset(dir_test):
    path_pdb = cache_5khb(dir_test)
    lines = run_unify_numbering(path_pdb, reset_initial_resid=False)
    line_hetatm = lines[140]
    line_hetatm_ref = (
        "HETATM    6  CG  FME A   1       0.765  -1.863  -1.806  1.00 35.55           C"
    )
    assert line_hetatm.strip() == line_hetatm_ref
    line_atom = lines[155]
    line_atom_ref = (
        "ATOM     21  CA  GLY A   2       2.192   3.731  -2.537  1.00 30.10           C"
    )
    assert line_atom.strip() == line_atom_ref


def test_unify_5khb_reset(dir_test):
    path_pdb = cache_5khb(dir_test)
    lines = run_unify_numbering(path_pdb, reset_initial_resid=True)
    line_hetatm = lines[140]
    line_hetatm_ref = (
        "HETATM    6  CG  FME A   1       0.765  -1.863  -1.806  1.00 35.55           C"
    )
    assert line_hetatm.strip() == line_hetatm_ref
    line_atom = lines[155]
    line_atom_ref = (
        "ATOM     21  CA  GLY A   2       2.192   3.731  -2.537  1.00 30.10           C"
    )
    assert line_atom.strip() == line_atom_ref
