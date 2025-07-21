# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

import os

from simlify import SimlifyConfig
from simlify.simulation.slurm import run_slurm_prep


def test_slurm_prep(dir_test):
    dir_work = os.path.join(dir_test, "tmp/slurm")
    path_slurm_write = "job.sbatch"
    simlify_config = SimlifyConfig()
    run_slurm_prep(dir_work, path_slurm_write, simlify_config)
