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

from simlify.structure.pdb.names import replace_residue_names
from simlify.structure.pdb.utils import keep_lines


def cache_5khb(dir_test: str) -> str:
    pdb_id = "5khb"
    path_pdb = os.path.join(dir_test, f"tmp/{pdb_id}.pdb")
    download_pdb(pdb_id, path_save=path_pdb)
    return path_pdb


def test_resname_5khb(dir_test):
    path_pdb = cache_5khb(dir_test)
    with open(path_pdb, "r") as f:
        lines_pdb = f.readlines()

    assert len(lines_pdb) == 7056
    lines_pdb = keep_lines(lines_pdb)
    assert len(lines_pdb) == 6901

    lines_pdb = lines_pdb[1:30]

    assert len(lines_pdb) == 29
    lines_pdb = replace_residue_names(lines_pdb, "FME", "A")
    assert len(lines_pdb) == 29

    line_ref = (
        "HETATM    7  SD  A   A   1      -0.649  -2.555  -2.729  1.00 33.41           S"
    )
    assert lines_pdb[6].strip() == line_ref
    line_ref = (
        "ATOM     22  C   GLY A   2       2.460   3.984  -4.007  1.00 50.22           C"
    )
    assert lines_pdb[21].strip() == line_ref
