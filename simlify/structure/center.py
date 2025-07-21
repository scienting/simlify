# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.


"""Module for centering molecular structures."""

from typing import Any

import MDAnalysis as mda
from loguru import logger

from simlify.structure.io import load_mda, write_mda
from simlify.structure.pos import get_com


def center_structure(u: mda.Universe) -> mda.Universe:
    """Centers an MDAnalysis Universe object by shifting its center of mass to the
    origin.

    Calculates the center of mass of all atoms in the provided MDAnalysis
    Universe and then subtracts this center of mass from the coordinates of
    each atom, effectively placing the system's center of mass at the coordinates
    (0.0, 0.0, 0.0).

    Args:
        u: The MDAnalysis Universe object representing the molecular system
            to be centered.

    Returns:
        The same MDAnalysis Universe object, but with the coordinates of all
        atoms shifted so that the center of mass is at the origin.
    """
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
    r"""Loads a molecular structure, centers it by its center of mass, and optionally
    saves the centered structure to a new file.

    This function orchestrates the process of reading a molecular structure
    using MDAnalysis, centering it using the
    [`center_structure`][simlify.structure.center.center_structure]
    function, and then writing the centered coordinates to a specified output file.

    Args:
        path_topo: The path to the topology file of the molecular system.
            Common file formats include PDB, GRO, and TPR.
        path_coords: The path to the coordinate file(s) of the molecular system.
            This can be a single file path (e.g., a trajectory file like TRR or XTC),
            a list of file paths, or `None` if the topology file already contains
            coordinate information (as is the case for some PDB files).
        path_output: The path to the output file where the centered structure
            will be saved. If this argument is `None`, the centered structure
            will not be written to a file. The file format is typically inferred
            from the file extension (e.g., `.pdb`, `.gro`).
        kwargs_universe: A dictionary of keyword arguments that will be passed
            to the [`load_mda`][simlify.structure.io.load_mda] function when
            creating the MDAnalysis Universe object. This allows for customization
            of how the structure is loaded (e.g., specifying atom names or residue
            numbers).
        kwargs_writer: A dictionary of keyword arguments that will be passed
            to the [`write_mda`][simlify.structure.io.write_mda] function when
            saving the centered structure. This can be used to control the output
            format or other writing options.
        overwrite: A boolean flag indicating whether to overwrite the output file
            if it already exists. If `True`, the existing file will be overwritten.
            If `False` and the file exists, an error might be raised depending
            on the `write_mda` function's behavior.

    Returns:
        The MDAnalysis Universe object representing the molecular system after
        it has been centered by its center of mass.

    Examples:
        Centering a PDB file and saving the result:

        ```python
        from simlify.structure.center import run_center_structure

        centered_universe = run_center_structure(
            path_topo="input.pdb",
            path_output="centered.pdb",
            overwrite=True,
        )
        print(f"Centered universe has {centered_universe.atoms.n_atoms} atoms.")
        ```

        Centering a topology and trajectory file without saving:

        ```python
        from simlify.structure.center import run_center_structure

        centered_universe = run_center_structure(
            path_topo="topol.gro",
            path_coords="traj.xtc",
        )
        print(
            f"Center of mass of the centered universe: {centered_universe.atoms.center_of_mass()}"
        )
        ```

        Passing additional arguments to the MDAnalysis Universe constructor:

        ```python
        from simlify.structure.center import run_center_structure

        centered_universe = run_center_structure(
            path_topo="protein.pdb",
            kwargs_universe={"atom_name_column": "name"},
            path_output="centered_protein.pdb",
        )
        ```
    """
    u = load_mda(path_topo, path_coords, **kwargs_universe)
    u = center_structure(u)

    if path_output:
        write_mda(u, u.atoms, path_output, overwrite=overwrite, **kwargs_writer)
    return u
