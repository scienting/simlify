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

from simlify.configs.slurm import SlurmConfig


def test_render_slurm():
    slurm_config = SlurmConfig()
    write_path = os.path.join(TMP_DIR, "job.slurm")
    slurm_config.write_render(write_path)
