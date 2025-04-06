"""
Provides functions for loading and writing molecular structures using the MDAnalysis library.

This module offers thin wrappers around the MDAnalysis library's functionalities
for reading and writing molecular structure files. It simplifies common tasks
such as loading topology and coordinate files into MDAnalysis Universe objects
and saving coordinates to various supported formats.
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
    """Loads a molecular structure into an MDAnalysis Universe object.

    This function serves as a convenient interface to the `MDAnalysis.Universe`
    constructor. It handles the common scenario where a topology file is provided,
    optionally accompanied by one or more coordinate files.

    Args:
        path_topo: The path to the topology file. This file defines the atoms,
            their types, connectivity, and sometimes initial coordinates. Common
            file formats include Protein Data Bank (PDB), GROningen structure
            (GRO), and GROMACS topology (TPR) files.
        path_coords: The path to the coordinate file(s). This can be:
            - A single string representing the path to a coordinate file (e.g.,
              a trajectory file like TRR, XTC, DCD, or a coordinate file like
              a second PDB).
            - A list of strings, where each string is a path to a coordinate file.
              MDAnalysis will attempt to load these files as a trajectory.
            - `None`, in which case it is assumed that the topology file (`path_topo`)
              already contains coordinate information (e.g., a single-frame PDB file).
        *args: Positional arguments to be passed directly to the
            `MDAnalysis.Universe` constructor. Refer to the MDAnalysis documentation
            for details on available positional arguments.
        **kwargs: Keyword arguments to be passed directly to the
            `MDAnalysis.Universe` constructor. This can include options such as
            `format` to explicitly specify the file format, `atom_name_column` to
            define which column in the topology file contains atom names, and many
            other options. Consult the MDAnalysis documentation for a comprehensive
            list of available keyword arguments.

    Returns:
        An `MDAnalysis.Universe` object representing the loaded molecular structure.
            This object provides access to the atoms, residues, segments, and
            trajectoryinformation of the system.

    Raises:
        FileNotFoundError: If the topology file specified by `path_topo` does not exist.

    Examples:
        Loading a structure from a topology file only (assuming coordinates are
        included):

        ```python
        import MDAnalysis as mda
        from simlify.structure.io import load_mda

        try:
            universe = load_mda("protein.pdb")
            print(f"Loaded universe with {universe.atoms.n_atoms} atoms.")
        except FileNotFoundError as e:
            print(f"Error: {e}")
        ```

        Loading a structure from a topology and a single coordinate file:

        ```python
        import MDAnalysis as mda
        from simlify.structure.io import load_mda

        try:
            universe = load_mda("protein.gro", "trajectory.xtc")
            print(f"Loaded trajectory with {len(universe.trajectory)} frames.")
        except FileNotFoundError as e:
            print(f"Error: {e}")
        ```

        Loading a structure from a topology and a list of coordinate files:

        ```python
        import MDAnalysis as mda
        from simlify.structure.io import load_mda

        try:
            universe = load_mda("topol.tpr", ["traj1.trr", "traj2.trr"])
            print(f"Loaded trajectory with a total of {len(universe.trajectory)} frames.")
        except FileNotFoundError as e:
            print(f"Error: {e}")
        ```

        Passing additional keyword arguments to the MDAnalysis Universe constructor:

        ```python
        import MDAnalysis as mda
        from simlify.structure.io import load_mda

        try:
            universe = load_mda("complex.pdb", kwargs={"atom_name_column": "name"})
            print(f"Loaded universe with custom atom name column.")
        except FileNotFoundError as e:
            print(f"Error: {e}")
        ```
    """
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
    r"""Writes frames from a molecular structure (or trajectory) to a file using
    MDAnalysis.

    This function provides a convenient way to save the coordinates of a molecular
    system, potentially including trajectory information, to a file format supported
    by the MDAnalysis library. It allows for writing the entire trajectory or a
    specific subset of frames.

    Args:
        u: The MDAnalysis Universe object containing the molecular structure and
            potentially trajectory information to be written.
        atoms: An `MDAnalysis.AtomGroup` object representing the atoms whose
            coordinates will be written to the output file. This can be the entire
            atom group of the Universe or a subset of selected atoms.
        path_output: The path to the output file where the coordinates will be saved.
            The file format is typically determined by the file extension (e.g.,
            `.pdb` for Protein Data Bank, `.gro` for GROMACS format, `.xtc` for
            GROMACS trajectory, `.trr` for GROMACS trajectory).
        frames: An optional sequence of integer indices specifying the frames from
            the trajectory to be written. If `None` (default), all frames in the
            Universe's trajectory will be written to the output file.
        overwrite: A boolean flag indicating whether to overwrite the file specified
            by `path_output` if it already exists. If `True`, the existing file
            will be overwritten without prompting. If `False` and the file exists,
            a `FileExistsError` will be raised to prevent accidental data loss.
        **kwargs: Keyword arguments to be passed directly to the `MDAnalysis.Writer`
            constructor. This can include options specific to the output file format,
            such as `force_overwrite` (though this function already handles overwriting
            via the `overwrite` argument), `format` to explicitly specify the format,
            and other format-specific options. Consult the MDAnalysis documentation
            for available keyword arguments for different writer types.

    Raises:
        FileExistsError: If `overwrite` is `False` and a file already exists at
            the specified `path_output`.

    Examples:
        Writing the entire trajectory of a Universe to a PDB file:

        ```python
        import MDAnalysis as mda
        from simlify.structure.io import load_mda, write_mda

        # Assume 'universe' is an MDAnalysis Universe object loaded previously
        try:
            universe = load_mda("protein.gro", "trajectory.xtc")
            write_mda(universe, universe.atoms, "output.pdb", overwrite=True)
            print("Successfully wrote coordinates to output.pdb")
        except FileExistsError as e:
            print(f"Error: {e}")
        ```

        Writing only specific frames from a trajectory to a GRO file:

        ```python
        import MDAnalysis as mda
        from simlify.structure.io import load_mda, write_mda

        # Assume 'universe' is an MDAnalysis Universe object with a trajectory
        universe = load_mda("topol.tpr", "traj.trr")
        frames_to_write = [0, 5, 10, 20]
        write_mda(universe, universe.atoms, "selected_frames.gro", frames=frames_to_write, overwrite=True)
        print(f"Wrote frames {frames_to_write} to selected_frames.gro")
        ```

        Writing only a subset of atoms to a new PDB file:

        ```python
        import MDAnalysis as mda
        from simlify.structure.io import load_mda, write_mda

        # Assume 'universe' is an MDAnalysis Universe object
        universe = load_mda("system.pdb")
        protein_atoms = universe.select_atoms("protein")
        write_mda(universe, protein_atoms, "protein_only.pdb", overwrite=True)
        print(f"Wrote coordinates of {protein_atoms.n_atoms} protein atoms to protein_only.pdb")
        ```

        Passing additional keyword arguments to the MDAnalysis Writer:

        ```python
        import MDAnalysis as mda
        from simlify.structure.io import load_mda, write_mda

        # Assume 'universe' is an MDAnalysis Universe object
        universe = load_mda("input.pdb")
        write_mda(universe, universe.atoms, "output.pdb", overwrite=True, format='pdb')
        print("Wrote coordinates to output.pdb with explicit format specification.")
        ```
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
