import os

from simlify import SimlifyConfig
from simlify.simulation.slurm import run_slurm_prep


def test_slurm_prep(dir_test):
    dir_work = os.path.join(dir_test, "tmp/slurm")
    path_slurm_write = "job.sbatch"
    simlify_config = SimlifyConfig()
    run_slurm_prep(dir_work, path_slurm_write, simlify_config)
