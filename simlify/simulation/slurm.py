import os

from simlify import SimlifyConfig


def run_slurm_prep(
    dir_work: str,
    path_slurm_write: str,
    simlify_config: SimlifyConfig,
) -> None:
    r"""Prepare files for simulations using slurm.

    Args:
        dir_work: Path to local directory where we will write the simulation files.
        path_slurm_write: Local path to write a slurm submission script with respect to
            `dir_write`.
        simlify_config: Simlify configuration.
    """
    os.makedirs(dir_work, exist_ok=True)

    lines_slurm = simlify_config.slurm.render(with_newlines=True)
    path_slurm_write = os.path.join(dir_work, path_slurm_write)
    with open(path_slurm_write, "w+", encoding="utf-8") as f:
        f.writelines(lines_slurm)
