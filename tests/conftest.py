# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.


import os
import urllib.request

import pytest

from simlify import SimlifyConfig, enable_logging
from simlify.schemas.amber import Amber22Schema
from simlify.structure.io import load_mda

TEST_DIR = os.path.dirname(__file__)
TMP_DIR = os.path.join(TEST_DIR, "tmp")
FILE_DIR = os.path.join(TEST_DIR, "files")


@pytest.fixture
def dir_test():
    return os.path.abspath(TEST_DIR)


@pytest.fixture(scope="session", autouse=True)
def turn_on_logging():
    enable_logging(10)


@pytest.fixture
def dir_structures():
    return os.path.join(TEST_DIR, "files/structures")


@pytest.fixture
def path_1jc0():
    return os.path.join(TEST_DIR, "files/structures/1JC0.pdb")


@pytest.fixture
def path_1jc0_prepped():
    return os.path.join(TEST_DIR, "files/structures/1JC0-prepped.pdb")


@pytest.fixture
def path_cro_fcrmod():
    return os.path.join(
        TEST_DIR, "files/ff/amber-ff-chromo-params/frcmod.xFPchromophores.2022"
    )


@pytest.fixture
def path_cro_lib():
    return os.path.join(
        TEST_DIR, "files/ff/amber-ff-chromo-params/xFPchromophores.lib.2022"
    )


@pytest.fixture
def amber_sim_standard_config():
    simlify_config = SimlifyConfig()
    simlify_config.engine = Amber22Schema()
    simlify_config.label = "01_min"
    simlify_config.run.dir_work = os.path.join(TEST_DIR, "tmp")
    simlify_config.engine.cli.compute_platform = "pmemd.MPI"
    simlify_config.engine.ff.protein = "ff19SB"
    simlify_config.engine.ff.water = "opc3"
    simlify_config.slurm.ntasks_per_node = 8
    simlify_config.engine.inputs.update(
        {
            "irest": 1,
            "ntx": 5,
            "ig": -1,
            "dt": 0.002,
            "nstlim": 500000,
            "nscm": 500,
            "ntr": 1,
            "restraint_wt": 0.5,
            "restraintmask": "!(:WAT) & (@C,CA,N,O,O5',P,O3',C3',C4',C5')",
            "ntb": 2,
            "ntf": 2,
            "ntc": 2,
            "cut": 10.0,
            "ntt": 3,
            "temp0": 300.0,
            "gamma_ln": 5.0,
            "ntp": 1,
            "barostat": 2,
            "pres0": 1.01325,
            "mcbarint": 100,
            "comp": 44.6,
            "taup": 1.0,
            "ntxo": 2,
            "ntwr": 5000,
            "ntpr": 500,
            "ntwx": 5000,
            "ioutfm": 1,
            "iwrap": 1,
        }
    )
    return simlify_config


def download_pdb(pdb_id, path_save=None):
    """Download and cache PDB files."""
    pdb_id = pdb_id.lower()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    path_save = path_save or os.path.join(TEST_DIR, f"tmp/{pdb_id}.pdb")

    if not os.path.exists(path_save):
        try:
            response = urllib.request.urlretrieve(url, path_save)
            print(response)
        except Exception as e:
            print(f"Error downloading PDB file: {e}")


@pytest.fixture
def u_1haj(dir_test):
    """Fixture to download and load the 1HAJ structure once per test module."""
    pdb_id = "1haj"
    path_pdb = os.path.join(dir_test, "tmp", f"{pdb_id}.pdb")
    download_pdb(pdb_id, path_save=path_pdb)
    return load_mda(path_pdb)
