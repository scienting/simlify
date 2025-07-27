# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

import os

from conftest import TMP_DIR

from simlify.schemas.slurm import SlurmSchema


def test_render_slurm():
    slurm_schema = SlurmSchema()
    write_path = os.path.join(TMP_DIR, "job.slurm")
    slurm_schema.write_render(write_path)
