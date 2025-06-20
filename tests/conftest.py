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


@pytest.fixture
def multichain_toy_pdb(dir_test):
    """Fixture to write a toy multichain PDB structure per test module"""
    pdb_contents = """HEADER    TEST MULTI-CHAIN PDB
ATOM      1  N   MET A   1      11.104  13.207  10.123  1.00 20.00           N  
ATOM      2  CA  MET A   1      12.560  13.500  10.123  1.00 20.00           C  
ATOM      3  C   MET A   1      13.104  12.750  11.333  1.00 20.00           C  
ATOM      4  O   MET A   1      12.500  11.780  11.678  1.00 20.00           O  
ATOM      5  N   GLY A   2      14.333  13.050  11.922  1.00 20.00           N  
ATOM      6  CA  GLY A   2      14.891  12.350  13.080  1.00 20.00           C  
ATOM      7  C   GLY A   2      14.011  12.700  14.317  1.00 20.00           C  
ATOM      8  O   GLY A   2      13.933  13.867  14.649  1.00 20.00           O  
TER
ATOM      9  N   ALA B   1      15.104  12.207  15.123  1.00 20.00           N  
ATOM     10  CA  ALA B   1      16.560  12.500  15.123  1.00 20.00           C  
ATOM     11  C   ALA B   1      17.104  11.750  16.333  1.00 20.00           C  
ATOM     12  O   ALA B   1      16.500  10.780  16.678  1.00 20.00           O  
ATOM     13  N   SER B   2      18.333  12.050  16.922  1.00 20.00           N  
ATOM     14  CA  SER B   2      18.891  11.350  18.080  1.00 20.00           C  
ATOM     15  C   SER B   2      18.011  11.700  19.317  1.00 20.00           C  
ATOM     16  O   SER B   2      17.933  12.867  19.649  1.00 20.00           O  
TER
END
"""
    path_pdb = os.path.join(dir_test, "tmp", "multichain_toy.pdb")
    with open(path_pdb, "w+") as f:
        f.write(pdb_contents)
    return path_pdb


@pytest.fixture
def multichain_hetatm_toy_pdb(dir_test):
    """Fixture to write a toy multichain PDB structure with 'HETATM' records per test module"""
    pdb_contents = """HEADER    TEST MULTI-CHAIN PDB
HETATM    1  N   MET A   1      11.104  13.207  10.123  1.00 20.00           N  
HETATM    2  CA  MET A   1      12.560  13.500  10.123  1.00 20.00           C  
ATOM      3  C   MET A   1      13.104  12.750  11.333  1.00 20.00           C  
ATOM      4  O   MET A   1      12.500  11.780  11.678  1.00 20.00           O  
ATOM      5  N   GLY A   2      14.333  13.050  11.922  1.00 20.00           N  
ATOM      6  CA  GLY A   2      14.891  12.350  13.080  1.00 20.00           C  
ATOM      7  C   GLY A   2      14.011  12.700  14.317  1.00 20.00           C  
ATOM      8  O   GLY A   2      13.933  13.867  14.649  1.00 20.00           O  
TER
HETATM    9  N   ALA B   1      15.104  12.207  15.123  1.00 20.00           N  
HETATM   10  CA  ALA B   1      16.560  12.500  15.123  1.00 20.00           C  
ATOM     11  C   ALA B   1      17.104  11.750  16.333  1.00 20.00           C  
ATOM     12  O   ALA B   1      16.500  10.780  16.678  1.00 20.00           O  
ATOM     13  N   SER B   2      18.333  12.050  16.922  1.00 20.00           N  
ATOM     14  CA  SER B   2      18.891  11.350  18.080  1.00 20.00           C  
ATOM     15  C   SER B   2      18.011  11.700  19.317  1.00 20.00           C  
ATOM     16  O   SER B   2      17.933  12.867  19.649  1.00 20.00           O  
TER
END
"""
    path_pdb = os.path.join(dir_test, "tmp", "multichain_hetatm_toy.pdb")
    with open(path_pdb, "w+") as f:
        f.write(pdb_contents)
    return path_pdb


@pytest.fixture
def multichain_seq_hetatm_toy_pdb(dir_test):
    """Fixture to write a toy multichain with sequential resid nunbering in the chains PDB structure per test module"""
    pdb_contents = """HEADER    TEST MULTI-CHAIN PDB
HETATM    1  N   MET A   1      11.104  13.207  10.123  1.00 20.00           N  
HETATM    2  CA  MET A   1      12.560  13.500  10.123  1.00 20.00           C  
ATOM      3  C   MET A   1      13.104  12.750  11.333  1.00 20.00           C  
ATOM      4  O   MET A   1      12.500  11.780  11.678  1.00 20.00           O  
ATOM      5  N   GLY A   2      14.333  13.050  11.922  1.00 20.00           N  
ATOM      6  CA  GLY A   2      14.891  12.350  13.080  1.00 20.00           C  
ATOM      7  C   GLY A   2      14.011  12.700  14.317  1.00 20.00           C  
ATOM      8  O   GLY A   2      13.933  13.867  14.649  1.00 20.00           O  
TER
HETATM    9  N   ALA B   3      15.104  12.207  15.123  1.00 20.00           N  
HETATM   10  CA  ALA B   3      16.560  12.500  15.123  1.00 20.00           C  
ATOM     11  C   ALA B   3      17.104  11.750  16.333  1.00 20.00           C  
ATOM     12  O   ALA B   3      16.500  10.780  16.678  1.00 20.00           O  
ATOM     13  N   SER B   4      18.333  12.050  16.922  1.00 20.00           N  
ATOM     14  CA  SER B   4      18.891  11.350  18.080  1.00 20.00           C  
ATOM     15  C   SER B   4      18.011  11.700  19.317  1.00 20.00           C  
ATOM     16  O   SER B   4      17.933  12.867  19.649  1.00 20.00           O  
TER
END
"""
    path_pdb = os.path.join(dir_test, "tmp", "multichain_seq_hetatm_toy.pdb")
    with open(path_pdb, "w+") as f:
        f.write(pdb_contents)
    return path_pdb


@pytest.fixture
def broken_multichain_toy_pdb(dir_test):
    """Fixture to write a toy multichain PDB structure where there is no 'TER' record and the numbering is nonsequential"""
    pdb_contents = """HEADER    TEST MULTI-CHAIN PDB
HETATM    1  N   MET A   1      11.104  13.207  10.123  1.00 20.00           N  
HETATM    2  CA  MET A   1      12.560  13.500  10.123  1.00 20.00           C  
ATOM      3  C   MET A   1      13.104  12.750  11.333  1.00 20.00           C  
ATOM      4  O   MET A   1      12.500  11.780  11.678  1.00 20.00           O  
ATOM      5  N   GLY A   2      14.333  13.050  11.922  1.00 20.00           N  
ATOM      6  CA  GLY A   2      14.891  12.350  13.080  1.00 20.00           C  
ATOM      7  C   GLY A   2      14.011  12.700  14.317  1.00 20.00           C  
ATOM      8  O   GLY A   2      13.933  13.867  14.649  1.00 20.00           O  
HETATM    9  N   ALA B   5      15.104  12.207  15.123  1.00 20.00           N  
HETATM   10  CA  ALA B   5      16.560  12.500  15.123  1.00 20.00           C  
ATOM     11  C   ALA B   5      17.104  11.750  16.333  1.00 20.00           C  
ATOM     12  O   ALA B   5      16.500  10.780  16.678  1.00 20.00           O  
ATOM     13  N   SER B   6      18.333  12.050  16.922  1.00 20.00           N  
ATOM     14  CA  SER B   6      18.891  11.350  18.080  1.00 20.00           C  
ATOM     15  C   SER B   6      18.011  11.700  19.317  1.00 20.00           C  
ATOM     16  O   SER B   6      17.933  12.867  19.649  1.00 20.00           O  
TER
END
"""
    path_pdb = os.path.join(dir_test, "tmp", "broken_multichain_hetatm_toy.pdb")
    with open(path_pdb, "w+") as f:
        f.write(pdb_contents)
    return path_pdb
