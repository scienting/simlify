# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

"""Module for optimizing the orientation of a molecular system to minimize its bounding box volume.

This module provides functions to rotate a set of atomic coordinates in 3D space
and to find the optimal rotation that results in the smallest possible bounding box
volume. It leverages the MDAnalysis library for reading molecular structures,
NumPy for numerical operations, SciPy for optimization and rotations, and loguru
for logging.

The core functionality includes:
- Rotating a set of atomic positions by specified Euler angles.
- Calculating the volume of the bounding box enclosing a set of positions.
- Using optimization algorithms to find the rotation angles that minimize this volume.
- A convenience function to load a PDB file, perform the optimization, and optionally
  save the rotated coordinates to a new PDB file.

This module is particularly useful for preparing molecular structures for simulations
or analyses where a compact representation or a specific orientation is desired.
"""

import MDAnalysis as mda
import numpy as np
import numpy.typing as npt
from loguru import logger
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation

from simlify.structure.dims import get_box_volume


def rotate_positions(
    positions: npt.NDArray[np.float64], rotation_v: npt.NDArray[np.float64]
) -> npt.NDArray[np.float64]:
    r"""Rotates a set of 3D Cartesian coordinates by applying Euler angles.

    This function takes an array of atomic positions and a vector of Euler angles
    (in degrees) and applies the corresponding rotation to the positions. The
    rotation is performed using the 'xyz' convention for Euler angles, meaning
    rotations are applied sequentially around the x, y, and z axes.

    Args:
        positions (npt.NDArray[np.float64]): A NumPy array of shape (N, 3) where N is the
            number of atoms, and each row represents the x, y, and z coordinates
            of an atom.
        rotation_v (npt.NDArray[np.float64]): A NumPy array of shape (3,) containing the
            Euler angles in degrees for rotation around the x, y, and z axes, respectively.

    Returns:
        npt.NDArray[np.float64]: A NumPy array of shape (N, 3) containing the rotated
            atomic coordinates.

    Examples:
        >>> import numpy as np
        >>> positions = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        >>> rotation_angles = np.array(
        ...     [90.0, 0.0, 0.0]
        ... )  # Rotate 90 degrees around the x-axis
        >>> rotated_positions = rotate_positions(positions, rotation_angles)
        >>> print(rotated_positions)
        [[ 1.00000000e+00  0.00000000e+00  0.00000000e+00]
         [ 0.00000000e+00  6.12323400e-17  1.00000000e+00]]
    """
    logger.trace("Applying this rotation: {}", rotation_v)
    rotation: Rotation = Rotation.from_euler("xyz", rotation_v, degrees=True)
    rotated_positions: npt.NDArray[np.float64] = rotation.apply(positions)
    return rotated_positions


def volume_objective_f(
    rotation_v: npt.NDArray[np.float64], positions: npt.NDArray[np.float64]
) -> float:
    r"""Objective function to be minimized, returning the volume of the bounding box.

    This function takes a set of Euler angles and a set of atomic positions as input.
    It first rotates the positions using the provided angles and then calculates the
    volume of the smallest axis-aligned bounding box that encloses the rotated positions.
    This function is designed to be used with optimization algorithms to find the
    rotation that minimizes the bounding box volume.

    Args:
        rotation_v (npt.NDArray[np.float64]): A NumPy array of shape (3,) containing the
            Euler angles in degrees for rotation around the x, y, and z axes, respectively.
        positions (npt.NDArray[np.float64]): A NumPy array of shape (N, 3) where N is the
            number of atoms, and each row represents the x, y, and z coordinates
            of an atom.

    Returns:
        float: The volume of the bounding box enclosing the rotated atomic positions.

    Notes:
        The bounding box volume is calculated using the `get_box_volume` function
        from the `simlify.structure.dims` module.

    Examples:
        >>> import numpy as np
        >>> positions = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        >>> rotation_angles = np.array([0.0, 0.0, 0.0])
        >>> volume = volume_objective_f(rotation_angles, positions)
        >>> print(f"{volume=}")
        volume=0.0
    """
    rotated_positions = rotate_positions(positions, rotation_v)
    return float(get_box_volume(rotated_positions))


def minimize_box(positions: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    r"""Rotates a system of particles to minimize the volume of its bounding box.

    This function takes the Cartesian coordinates of a system of particles and uses
    numerical optimization to find the Euler angles that, when applied to the system,
    result in the smallest possible axis-aligned bounding box volume. The optimization
    is performed using the `scipy.optimize.minimize` function with the
    `volume_objective_f` as the objective function. The optimization is bounded to
    rotation angles between 0 and 360 degrees for each axis.

    Args:
        positions (npt.NDArray[np.float64]): A NumPy array of shape (N, 3) where N is the
            number of particles, and each row represents the x, y, and z coordinates.

    Returns:
        npt.NDArray[np.float64]: A NumPy array of shape (3,) containing the optimized
            Euler angles (in degrees) that result in the minimal bounding box volume.

    Examples:
        >>> import numpy as np
        >>> positions = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        >>> optimized_angles = minimize_box(positions)
        >>> print(f"{optimized_angles=}")
        optimized_angles=array([0., 0., 0.])
    """
    logger.info("Minimizing box size")
    initial_angles = np.zeros((3,), dtype=np.float64)
    logger.trace("Initial angles: {}", initial_angles)
    results = minimize(
        volume_objective_f,
        initial_angles,
        args=(positions,),
        bounds=[(0, 360)] * 3,
    )
    rotated_positions: npt.NDArray[np.float64] = results.x
    return rotated_positions


def run_minimize_box(
    pdb_path: str, output_path: str | None = None
) -> npt.NDArray[np.float64]:
    r"""Minimizes the bounding box size of a molecular system described in a PDB file by rotating it.

    This function serves as a high-level interface to load a molecular structure from a PDB
    file, optimize its orientation to minimize the bounding box volume using the
    `minimize_box` function, and optionally save the rotated coordinates to a new PDB file.

    Args:
        pdb_path (str): The path to the input PDB file containing the molecular structure.
        output_path (str | None, optional): The path to save a new PDB file with the rotated
            coordinates. If `None`, no file is written, and only the optimized atomic
            positions are returned. Defaults to `None`.

    Returns:
        npt.NDArray[np.float64]: A NumPy array of shape (N, 3) containing the atomic
            positions that have been rotated to achieve the minimum bounding box volume.

    Raises:
        FileNotFoundError: If the specified `pdb_path` does not exist.
        IOError: If there is an error reading the PDB file or writing to the output file.

    Examples:
        To load a PDB file named "input.pdb", minimize its bounding box volume, and save the result to "optimized.pdb":

        >>> optimized_positions = run_minimize_box(
        ...     "input.pdb", output_path="optimized.pdb"
        ... )

        To perform the minimization and get the optimized positions without saving to a file:

        >>> optimized_positions = run_minimize_box("input.pdb")
        >>> print(optimized_positions)
        [[...], [...], ...]
    """
    u = mda.Universe(pdb_path)
    optimized_rotation_v = minimize_box(u.atoms.positions)
    optimized_positions = rotate_positions(u.atoms.positions, optimized_rotation_v)

    if isinstance(output_path, str):
        u.atoms.positions = optimized_positions
        u.atoms.write(output_path)
    return optimized_positions
