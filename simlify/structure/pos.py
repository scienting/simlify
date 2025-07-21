# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

"""
Provides functions for interacting with the positions of atoms in molecular structures.

This module offers utilities for analyzing and manipulating the spatial coordinates
of atoms within an MDAnalysis Universe. Currently, it includes functionality to
calculate the center of mass of a selected group of atoms.
"""

import MDAnalysis as mda
import numpy as np
import numpy.typing as npt
from loguru import logger


def get_com(atoms: mda.AtomGroup) -> npt.NDArray[np.float64]:
    r"""Compute the center of mass for a group of atoms.

    The center of mass (COM) is the average position of all the parts of the
    system, weighted by their masses. For a collection of atoms, it is calculated
    as the weighted average of their coordinates, where the weights are the masses
    of the atoms. This function leverages the built-in center of mass calculation
    provided by MDAnalysis for an `AtomGroup`.

    Args:
        atoms: An MDAnalysis `AtomGroup` object representing the selection of
            atoms for which the center of mass needs to be computed. This object
            contains information about the atoms, including their positions and
            masses (if available in the topology).

    Returns:
        A NumPy array of shape (3,) containing the x, y, and z coordinates of
        the center of mass of the provided `AtomGroup`. The units of these
        coordinates will be the same as the units of the atomic coordinates
        within the MDAnalysis Universe (typically Angstroms).

    Examples:
        Calculating the center of mass of all atoms in a loaded structure:

        ```python
        import MDAnalysis as mda
        from simlify.structure.io import load_mda
        from simlify.structure.pos import get_com

        # Load a molecular structure
        try:
            universe = load_mda("protein.pdb")
            all_atoms = universe.atoms

            # Calculate the center of mass of all atoms
            center_of_mass = get_com(all_atoms)
            print(f"Center of mass of all atoms: {center_of_mass}")
        except FileNotFoundError as e:
            print(f"Error: {e}")
        ```

        Calculating the center of mass of a specific selection of atoms (e.g., protein):

        ```python
        import MDAnalysis as mda
        from simlify.structure.io import load_mda
        from simlify.structure.pos import get_com

        # Load a molecular structure
        try:
            universe = load_mda("complex.gro", "traj.xtc")
            protein = universe.select_atoms("protein")

            # Calculate the center of mass of the protein atoms
            protein_com = get_com(protein)
            print(f"Center of mass of the protein: {protein_com}")
        except FileNotFoundError as e:
            print(f"Error: {e}")
        ```
    """
    logger.info("Computing center of mass (com)")
    com: npt.NDArray[np.float64] = atoms.center_of_mass()
    logger.debug("center of mass: {}", com)
    return com
