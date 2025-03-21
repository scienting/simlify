from typing import Any

import os

import MDAnalysis as mda
from loguru import logger


def run_select_atoms(
    select: str,
    path_topo: str,
    path_output: str | None = None,
    path_coords: str | list[str] | None = None,
    kwargs_universe: dict[str, Any] = dict(),
    kwargs_writer: dict[str, Any] = dict(),
    overwrite: bool = False,
) -> mda.core.groups.AtomGroup:
    r"""Select atoms from molecular structure using
    [MDAnalysis](https://docs.mdanalysis.org/stable) and optionally write new
    coordinates.

    Args:
        select: [MDAnalysis selection string](https://docs.mdanalysis.org/stable/documentation_pages/selections.html).
        path_topo: Path to topology file.
        path_output: Path to save new coordinate file. If `None`, then no file is
            written.
        path_coords: Paths to coordinate file(s) to load.
        kwargs_universe: Keyword arguments for initializing the
            [MDAnalysis Universe](https://docs.mdanalysis.org/stable/documentation_pages/core/universe.html#MDAnalysis.core.universe.Universe).
        kwargs_writer: Keyword arguments for
            [MDAnalysis writers](https://docs.mdanalysis.org/stable/documentation_pages/coordinates/init.html#writers).
        overwrite: Overwrites coordinate file located at `path_output`. This will error
            out if `False` and the file exists.

    Returns:
        Atomic positions after atom selection.

    Example:
        Suppose we want to extract just protein atoms from an Amber simulation.

        ```python
        run_select_atoms(
            select="protein",
            path_topo="mol.prmtop",
            path_output="protein.nc",
            path_coords="traj.nc",
        )
        ```
    """
    if not os.path.exists(path_topo):
        logger.critical(f"Could not find topology file at {path_topo}")
        raise

    u = mda.Universe(path_topo, path_coords, **kwargs_universe)

    atoms = u.select_atoms(select)

    if path_output:
        if os.path.exists(path_output):
            if not overwrite:
                logger.critical(
                    f"Overwrite is False and file already exists at {path_output}"
                )

        with mda.Writer(path_output, atoms.n_atoms) as W:
            for ts in u.trajectory:
                W.write(atoms, **kwargs_writer)
    return atoms
