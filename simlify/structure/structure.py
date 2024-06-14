import MDAnalysis as mda
import numpy as np
import numpy.typing as npt
from loguru import logger


def get_box_lengths(positions: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """
    Compute the lengths of the box edges.

    This function calculates the lengths of the edges of a box that
    encompasses a set of atomic positions by finding the difference
    between the maximum and minimum coordinates in each dimension.

    Args:
        positions: A 2D array of shape (N, 3) containing the Cartesian coordinates of
            N atoms.

    Returns:
        A 1D array of shape (3,) containing the lengths
            of the box edges in the x, y, and z dimensions, respectively.

    Examples:
        ```python
        >>> positions = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]])
        >>> print(get_box_lengths(positions))
        array([2., 2., 2.])
        ```
    """
    logger.trace("Computing box lengths")
    box_lengths: npt.NDArray[np.float64] = np.max(positions, axis=0) - np.min(
        positions, axis=0
    )
    logger.trace("Box lengths: {}", box_lengths)
    return box_lengths


def get_box_vectors(positions: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """
    Compute the vectors of the box edges.

    This function calculates the vectors of the edges of a box that
    encompasses a set of atomic positions. The vectors are represented
    as a diagonal matrix where the diagonal elements correspond to the
    lengths of the box edges.

    Args:
        positions: A 2D array of shape (N, 3) containing the Cartesian coordinates of
            N atoms.

    Returns:
        A 2D array of shape (3, 3) where the diagonal elements represent the vectors
            of the box edges.

    Examples:
        ```python
        >>> positions = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]])
        >>> print(get_box_vectors(positions))
        array([[2., 0., 0.],
               [0., 2., 0.],
               [0., 0., 2.]])
        ```
    """
    box_lengths = get_box_lengths(positions)

    logger.trace("Computing box vectors")
    box_vectors = np.zeros((3, 3), dtype=np.float64)
    np.fill_diagonal(box_vectors, box_lengths)
    logger.trace("Box vector:\n{}", box_vectors)

    return box_vectors


def get_com(universe: mda.Universe) -> npt.NDArray[np.float64]:
    """
    Compute the center of mass of the universe.

    This function calculates the center of mass (COM) of all atoms in a given
    MDAnalysis universe object.

    Args:
        universe: An MDAnalysis Universe object containing
            the atomic positions and masses.

    Returns:
        A 1D array of shape (3,) containing the Cartesian coordinates
            of the center of mass.

    Examples:
        ```python
        >>> u = mda.Universe('topology.psf', 'trajectory.dcd')
        >>> print(get_com(u))
        array([10.0, 10.0, 10.0])
        ```
    """
    logger.info("Computing center of mass (com)")
    com: npt.NDArray[np.float64] = universe.atoms.center_of_mass()
    logger.debug("com: {}", com)
    return com


def get_box_volume(positions: npt.NDArray[np.float64]) -> float:
    """
    Compute the volume of the box.

    This function calculates the volume of a box that encompasses a set of atomic
    positions by computing the determinant of the matrix formed by the box vectors.

    Args:
        positions: A 2D array of shape (N, 3) containing
            the Cartesian coordinates of N atoms.

    Returns:
        The volume of the box.

    Examples:
        ```python
        >>> positions = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]])
        >>> print(get_box_volume(positions))
        8.0
        ```
    """
    box_vectors = get_box_vectors(positions)
    logger.trace("Computing box volume")
    volume: float = np.linalg.det(box_vectors)
    logger.trace("Box volume: {}", volume)
    return volume
