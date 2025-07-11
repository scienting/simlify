from typing import Any, Sequence

import MDAnalysis as mda
from loguru import logger

from simlify.structure.io import load_mda, write_mda


def extract_atoms(
    path_topo: str,
    select: str | None = None,
    frames: Sequence[int] | None = None,
    path_output: str | None = None,
    path_coords: str | list[str] | None = None,
    kwargs_universe: dict[str, Any] = dict(),
    kwargs_writer: dict[str, Any] = dict(),
    overwrite: bool = False,
) -> mda.core.groups.AtomGroup:
    r"""Extract atoms or frames from molecular structure(s) using
    [MDAnalysis](https://docs.mdanalysis.org/stable).

    Args:
        path_topo: Path to topology file.
        select: [MDAnalysis selection string](https://docs.mdanalysis.org/stable/documentation_pages/selections.html).
        frames: Extract frames from a trajectory. If `None`, the whole trajectory will
            be included. This is only relevant when `path_output` is specified.
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
            path_topo="mol.prmtop",
            select="protein",
            frames=[0, -1],
            path_output="protein.nc",
            path_coords="traj.nc",
        )
        ```
    """
    if (select is None) and (frames is None):
        logger.warning(
            "Both selection and frames are None. No changes to the structure(s) will be made."
        )

    u: mda.Universe = load_mda(path_topo, path_coords, **kwargs_universe)

    if select:
        logger.info(f"Making selection: {select}")
        atoms = u.select_atoms(select)
    else:
        logger.info("Keeping all atoms")
        atoms = u.atoms

    if path_output:
        write_mda(u, atoms, path_output, frames, overwrite)
    return atoms
