"""
Interacting with atom positions.
"""

import MDAnalysis as mda
import numpy as np
import numpy.typing as npt
from loguru import logger


def get_com(atoms: mda.AtomGroup) -> npt.NDArray[np.float64]:
    r"""Compute the center of mass.

    Args:
        atoms: Atoms object from MDAnalysis.

    Returns:
        Center of mass of atoms group.
    """
    logger.info("Computing center of mass (com)")
    com: npt.NDArray[np.float64] = atoms.center_of_mass()
    logger.debug("center of mass: {}", com)
    return com
