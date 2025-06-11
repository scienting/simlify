import os

from simlify.schemas.slurm import SlurmSchema

from conftest import TMP_DIR


def test_render_slurm():
    slurm_schema = SlurmSchema()
    write_path = os.path.join(TMP_DIR, "job.slurm")
    slurm_schema.write_render(write_path)
