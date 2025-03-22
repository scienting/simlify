from typing import Any, Sequence

import os

import MDAnalysis as mda
from loguru import logger


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
    if not os.path.exists(path_topo):
        logger.critical(f"Could not find topology file at {path_topo}")
        raise

    if (select is None) and (frames is None):
        logger.warning(
            "Both selection and frames are None. No changes to the structure(s) will be made."
        )

    if path_coords:
        logger.debug("Loading all coordinate files")
        u = mda.Universe(path_topo, path_coords, **kwargs_universe)
    else:
        logger.debug("No coordinate file specified.")
        u = mda.Universe(path_topo, **kwargs_universe)

    if select:
        logger.debug(f"Making selection: {select}")
        atoms = u.select_atoms(select)
    else:
        logger.debug("Keeping all atoms")
        atoms = u.atoms

    if path_output:
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
                    W.write(atoms, **kwargs_writer)
            else:
                logger.info("Writing all frames")
                for ts in u.trajectory:
                    W.write(atoms, **kwargs_writer)
    return atoms
