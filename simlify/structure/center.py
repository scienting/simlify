from typing import Any

import MDAnalysis as mda
from loguru import logger

from simlify.structure.io import load_mda, write_mda
from simlify.structure.pos import get_com


def center_structure(u: mda.Universe) -> mda.Universe:
    com = get_com(u.atoms)
    # Shift the entire system so that the COM is at the origin (0, 0, 0)
    logger.info("Centering structure")
    u.atoms.positions -= com
    return u


def run_center_structure(
    path_topo: str,
    path_coords: str | list[str] | None = None,
    path_output: str | None = None,
    kwargs_universe: dict[str, Any] = dict(),
    kwargs_writer: dict[str, Any] = dict(),
    overwrite: bool = False,
) -> mda.Universe:
    r"""Center structure by it's center of mass.

    Args:
        path_topo: Path to PDB file.
        output_path: Path to save new PDB file. If `None`, then no file is written.

    Returns:
        Universe after centering.
    """
    u = load_mda(path_topo, path_coords, **kwargs_universe)
    u = center_structure(u)

    if path_output:
        write_mda(u, u.atoms, path_output, overwrite=overwrite, **kwargs_writer)
    return u
