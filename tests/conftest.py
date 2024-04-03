import os
from datetime import timedelta

import pytest
from pyarrow import fs

from simlify import enable_logging
from simlify.simulation.amber.contexts import AMBER_PROTEIN_STANDARD_CONTEXT
from simlify.simulation.contexts import SimContextManager

TEST_DIR = os.path.dirname(__file__)
GCS_CACHE_DIR = os.path.join(TEST_DIR, "files/gcs-cache/")
os.makedirs(GCS_CACHE_DIR, exist_ok=True)

gcs = fs.GcsFileSystem(anonymous=True, retry_time_limit=timedelta(seconds=5))


@pytest.fixture
def test_dir():
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
def amber_protein_standard_context():
    context_manager = SimContextManager(**AMBER_PROTEIN_STANDARD_CONTEXT)
    return context_manager


@pytest.fixture
def amber_simulation_standard_context(amber_protein_standard_context):
    context_manager = amber_protein_standard_context
    context_manager.name_stage = "01_min"
    context_manager.dir_input = TEST_DIR
    context_manager.dir_output = TEST_DIR
    context_manager.path_topo = os.path.join(TEST_DIR, "mol.prmtop")
    context_manager.path_input = os.path.join(
        TEST_DIR, f"{context_manager.name_stage}.in"
    )
    context_manager.path_restart_prev = os.path.join(TEST_DIR, "mol.inpcrd")
    context_manager.path_coord_ref = os.path.join(TEST_DIR, "mol.inpcrd")
    context_manager.compute_platform = "mpi"
    context_manager.cpu_cores = 12
    context_manager.input_kwargs = {
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
    return context_manager


@pytest.fixture
def dir_amber_sims():
    return os.path.join(TEST_DIR, "files/simulations/amber")


@pytest.fixture
def gcs_fs():
    return gcs


@pytest.fixture
def gcs_cache_dir():
    return GCS_CACHE_DIR


@pytest.fixture
def gcs_uuid():
    return "f7498a8c-d021-491c-a343-10151e81434a"
