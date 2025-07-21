# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

"""Module for preparing files necessary to run simulations using the SLURM workload manager."""

import os

from simlify import SimlifyConfig


def run_slurm_prep(
    dir_work: str,
    path_slurm_write: str,
    simlify_config: SimlifyConfig,
) -> None:
    r"""Prepares the necessary files for running simulations using the SLURM workload
    manager.

    This function takes a working directory, a path for the SLURM submission script, and
    the Simlify configuration. It ensures that the working directory exists and then renders
    a SLURM submission script based on the SLURM configuration found within the provided
    `simlify_config`. This script is then written to the specified path within the working
    directory.

    Args:
        dir_work: The absolute or relative path to the local directory where the
            simulation files, including the SLURM submission script, will be written.
            This directory will be created if it does not already exist.
        path_slurm_write: The path to write the SLURM submission script, relative to
            the `dir_work`. For example, if `dir_work` is `/home/user/simulations` and
            `path_slurm_write` is `submit.sh`, the script will be written to
            `/home/user/simulations/submit.sh`.
        simlify_config: An instance of the Simlify configuration object,
            which should contain a `slurm` attribute with the configuration details for
            generating the SLURM submission script.

    Raises:
        OSError: If there is an error creating the working directory or writing the
            SLURM submission script.

    Examples:
        To prepare SLURM files for a simulation in the directory "my_simulation" and
        write the submission script to "submit.sh":

        >>> from simlify import SimlifyConfig
        >>> config = SimlifyConfig()
        >>> run_slurm_prep(
        ...     dir_work="my_simulation",
        ...     path_slurm_write="submit.sh",
        ...     simlify_config=config,
        ... )
    """
    os.makedirs(dir_work, exist_ok=True)

    lines_slurm = simlify_config.slurm.render(with_newlines=True)
    path_slurm_write = os.path.join(dir_work, path_slurm_write)
    with open(path_slurm_write, "w+", encoding="utf-8") as f:
        f.writelines(lines_slurm)
