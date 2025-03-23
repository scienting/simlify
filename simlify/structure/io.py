"""
Loading structures
"""

from typing import Any, Sequence

import os

import MDAnalysis as mda
from loguru import logger


def load_mda(
    path_topo: str,
    path_coords: str | list[str] | None = None,
    *args: tuple[Any, ...],
    **kwargs: dict[str, Any],
) -> mda.Universe:
    """Thin wrapper around MDAnalysis Universe."""
    if not os.path.exists(path_topo):
        logger.critical(f"Could not find topology file at {path_topo}")
        raise

    if path_coords:
        logger.debug("Loading all coordinate files")
        u = mda.Universe(path_topo, path_coords, *args, **kwargs)
    else:
        logger.debug("No coordinate file specified.")
        u = mda.Universe(path_topo, *args, **kwargs)

    return u


def write_mda(
    u: mda.Universe,
    atoms: mda.AtomGroup,
    path_output: str,
    frames: Sequence[int] | None = None,
    overwrite: bool = False,
    **kwargs: dict[str, Any],
) -> None:
    r"""Write frames from molecular structure(s) using
    [MDAnalysis](https://docs.mdanalysis.org/stable).

    Args:
        u: MDAnalysis universe.
        atoms: Selected (or all) atoms of the same `u` Universe.
        path_output: Path to save new coordinate file.
        frames: Save only specific frames from a trajectory. If `None`, the whole
            trajectory will be included.
        overwrite: Overwrites coordinate file located at `path_output`. This will error
            out if `False` and the file exists.
    """
    if os.path.exists(path_output):
        logger.info(f"File at {path_output} already exists")
        if not overwrite:
            logger.critical(
                f"Overwrite is False and file already exists at {path_output}"
            )
        else:
            logger.info(f"Will overwrite {path_output}")

    with mda.Writer(path_output, atoms.n_atoms) as W:
        logger.info(f"Writing coordinates to {path_output}")
        if frames:
            logger.info(f"Writing {len(frames)} frames")
            for frame in frames:
                u.trajectory[frame]
                W.write(atoms, **kwargs)
        else:
            logger.info("Writing all frames")
            for ts in u.trajectory:
                W.write(atoms, **kwargs)
