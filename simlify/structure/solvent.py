from loguru import logger

from ..simulation.contexts import SimContextManager


def get_ion_counts(
    sim_context_manager: SimContextManager,
    charge_net: int | float,
    n_waters: int,
    water_molecule_volume: float = 28.78277638661025,
) -> dict[str, int]:
    r"""Compute the number of cations and anions to achieve desired ionic strength.

    Args:
        sim_context_manager: A simulation context for system preparation.
        charge_net: Total system charge. Will add counter ions if `charge_neutralize`
            is `True` in `sim_context_manager`.
        n_waters: Total number of water molecules in the system.
        water_molecule_volume: Approximate volume of a water molecule in Å<sup>3</sup>.

    Returns:
        Cations and anions that needs to be added to the system. For example, `2` would
        mean two cations need to be added, but `-2` would mean two anions.
    """
    logger.info("Computing number of extra ions")
    context = sim_context_manager.get()
    water_box_volume: float = n_waters * water_molecule_volume  # A^3
    logger.debug("Volume of water: {} A^3", water_box_volume)
    water_box_volume /= 1e27  # L
    n_ions = context["solvent_ionic_strength"] * water_box_volume  # moles
    n_ions *= 6.0221409e23  # atoms
    extra_ions = int(round(n_ions, 0))
    ions = {
        "charge_cation_num": context["charge_cation_extra"],
        "charge_anion_num": context["charge_anion_extra"],
    }
    ions = {k: v + extra_ions for k, v in ions.items()}
    if context["charge_neutralize"]:
        if charge_net < 0:
            ions["charge_cation_num"] += abs(int(charge_net))
        elif charge_net > 0:
            ions["charge_anion_num"] += abs(int(charge_net))
    logger.debug("Ions to add {}", ions)
    return ions
