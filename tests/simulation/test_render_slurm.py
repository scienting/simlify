import os

from conftest import TMP_DIR

from simlify.configs.slurm import SlurmConfig


def test_render_slurm():
    slurm_config = SlurmConfig()
    write_path = os.path.join(TMP_DIR, "job.slurm")
    slurm_config.write_render(write_path)
