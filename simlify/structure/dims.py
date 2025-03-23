"""
Dimensions of systems.

This module provides utility functions to calculate geometric properties
of a molecular system such as the bounding box lengths, box vectors,
and the volume of the box required to enclose the atoms.
"""

import numpy as np
import numpy.typing as npt
from loguru import logger


def get_box_lengths(positions: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    r"""Compute lengths of box edges.

    This function calculates the lengths of the edges of a box that
    encompasses a set of atomic positions by finding the difference
    between the maximum and minimum coordinates in each dimension.

    Args:
        positions: A 2D NumPy array of shape (N, 3), where N is the
            number of atoms and each row represents the Cartesian coordinates
            [x, y, z] of an atom.

    Returns:
        A 1D NumPy array of shape (3,) containing the lengths of the
            box in the x, y, and z directions.

    Example:
        ```python
        >>> import numpy as np
        >>> from box_utils import get_box_lengths
        >>> positions = np.array([[0, 0, 0], [1, 2, 3]])
        >>> get_box_lengths(positions)
        array([1., 2., 3.])
        ```
    """
    logger.trace("Computing box lengths")
    box_lengths: npt.NDArray[np.float64] = np.max(positions, axis=0) - np.min(
        positions, axis=0
    )
    logger.trace("Box lengths: {}", box_lengths)
    return box_lengths


def get_box_vectors(positions: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    r"""Construct orthogonal box vectors.

    This function calculates the vectors of the edges of a box that
    encompasses a set of atomic positions. The vectors are represented
    as a diagonal matrix where the diagonal elements correspond to the
    lengths of the box edges.

    Args:
        positions: A 2D NumPy array of shape (N, 3), where N is the
            number of atoms and each row represents the Cartesian coordinates
            [x, y, z] of an atom.

    Returns:
        A 3x3 NumPy array representing the box vectors in matrix form.
            The matrix is diagonal since the box is assumed to be orthorhombic.

    Example:
        ```python
        >>> import numpy as np
        >>> from box_utils import get_box_vectors
        >>> positions = np.array([[0, 0, 0], [1, 2, 3]])
        >>> get_box_vectors(positions)
        array([[1., 0., 0.],
               [0., 2., 0.],
               [0., 0., 3.]])
        ```
    """
    box_lengths = get_box_lengths(positions)

    logger.trace("Computing box vectors")
    box_vectors = np.zeros((3, 3), dtype=np.float64)
    np.fill_diagonal(box_vectors, box_lengths)
    logger.trace("Box vector:\n{}", box_vectors)

    return box_vectors


def get_box_volume(
    positions: npt.NDArray[np.float64],
) -> float:
    r"""Calculate volume of the enclosing orthorhombic box.

    Computes the volume of the orthorhombic box that encapsulates all atoms.
    The box vectors are assumed to form an orthogonal basis.

    Args:
        positions: A 2D NumPy array of shape (N, 3), where N is the
            number of atoms and each row represents the Cartesian coordinates
            [x, y, z] of an atom.

    Returns:
        Volume of the enclosing box.

    Example:
        ```python
        >>> import numpy as np
        >>> from box_utils import get_box_volume
        >>> positions = np.array([[0, 0, 0], [1, 2, 3]])
        >>> get_box_volume(positions)
        6.0
        ```
    """
    box_vectors = get_box_vectors(positions)
    logger.trace("Computing box volume")
    volume = float(np.linalg.det(box_vectors))
    logger.trace("Box volume: {}", volume)
    return volume
