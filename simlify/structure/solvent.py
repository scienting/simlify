"""
Provides methods to interact or prepare solvents for molecular simulation.
"""

from loguru import logger

from simlify.configs import SimlifyConfig


def get_ion_counts(
    simlify_config: SimlifyConfig,
    charge_net: int | float,
    n_waters: int,
    water_molecule_volume: float = 28.78277638661025,
) -> dict[str, int]:
    r"""Computes the number of cations and anions to add to a molecular simulation
    system to achieve a desired ionic strength and/or neutralize the system's net
    charge.

    This function calculates the number of ions needed based on the target ionic
    strength specified in the `SimlifyConfig`, the estimated volume of the water
    solvent, and the net charge of the solute. It also considers any pre-defined
    extra ions specified in the `SimlifyConfig`.

    The volume of the solvent is estimated based on the number of water molecules
    and an approximate volume per water molecule. The function then calculates
    the number of ions required to reach the target ionic strength in this volume.

    If the `charge_neutralize` flag is set to `True` in the `SimlifyConfig`, the
    function will also add counterions to neutralize the net charge of the system
    (`charge_net`). If the net charge is negative, cations will be added; if positive,
    anions will be added.

    Args:
        simlify_config: An instance of the `SimlifyConfig` class containing
            simulation parameters, including the desired solvent ionic strength
            (`simlify_config.solution.solvent_ionic_strength`), the number of
            extra cations (`simlify_config.solution.charge_cation_extra`), the
            number of extra anions (`simlify_config.solution.charge_anion_extra`),
            and the flag for charge neutralization
            (`simlify_config.solution.charge_neutralize`).
        charge_net: The net charge of the solute (e.g., a protein or DNA molecule)
            in the system before any ions are added. This value is used to determine
            the number of counterions needed for neutralization if specified in the
            `simlify_config`.
        n_waters: The total number of water molecules intended to be present in the
            simulation system. This number is used to estimate the volume of the
            solvent.
        water_molecule_volume: An approximate volume of a single water molecule in
            Angstroms cubed (Å³). The default value is based on the density of
            water at standard conditions.

    Returns:
        A dictionary containing the calculated number of cations and anions to be
        added to the system. The dictionary has the following keys:
        - `"charge_cation_num"`: The total number of cations to add, including
          any pre-defined extra cations and those required for ionic strength
          and/or charge neutralization.
        - `"charge_anion_num"`: The total number of anions to add, including
          any pre-defined extra anions and those required for ionic strength
          and/or charge neutralization.

    Examples:
        Calculating ion counts for a system with a net positive charge and a target
        ionic strength:

        ```python
        from simlify.configs import SimlifyConfig
        from simlify.solvate.ions import get_ion_counts

        # Create a SimlifyConfig instance
        config = SimlifyConfig(
            solution=dict(
                solvent_ionic_strength=0.150,  # 150 mM ionic strength
                charge_cation_extra=0,
                charge_anion_extra=0,
                charge_neutralize=True,
            )
        )
        net_charge = 2  # Example net positive charge
        num_waters = 10000

        ion_counts = get_ion_counts(config, net_charge, num_waters)
        print(f"Number of ions to add: {ion_counts}")
        # Expected output might vary slightly due to rounding, but will include anions
        # for neutralization.
        ```

        Calculating ion counts for a system with a net negative charge and no charge
        neutralization:

        ```python
        from simlify.configs import SimlifyConfig
        from simlify.solvate.ions import get_ion_counts

        # Create a SimlifyConfig instance with charge neutralization disabled
        config = SimlifyConfig(
            solution=dict(
                solvent_ionic_strength=0.100,  # 100 mM ionic strength
                charge_cation_extra=1,
                charge_anion_extra=0,
                charge_neutralize=False,
            )
        )
        net_charge = -3  # Example net negative charge
        num_waters = 5000

        ion_counts = get_ion_counts(config, net_charge, num_waters)
        print(f"Number of ions to add: {ion_counts}")
        # Expected output will include the pre-defined extra cation but no ions for
        # neutralization.
        ```

        Calculating ion counts for a neutral system with a specific ionic strength:

        ```python
        from simlify.configs import SimlifyConfig
        from simlify.solvate.ions import get_ion_counts

        # Create a SimlifyConfig instance for a neutral system
        config = SimlifyConfig(
            solution=dict(
                solvent_ionic_strength=0.200,  # 200 mM ionic strength
                charge_cation_extra=0,
                charge_anion_extra=0,
                charge_neutralize=False,
            )
        )
        net_charge = 0
        num_waters = 15000

        ion_counts = get_ion_counts(config, net_charge, num_waters)
        print(f"Number of ions to add: {ion_counts}")
        # Expected output will show the number of cations and anions needed for the
        # target ionic strength.
        ```
    """
    logger.info("Computing number of extra ions")
    water_box_volume: float = n_waters * water_molecule_volume  # A^3
    logger.debug("Volume of water: {} A^3", water_box_volume)
    water_box_volume /= 1e27  # L
    n_ions = simlify_config.solution.solvent_ionic_strength * water_box_volume  # moles
    n_ions *= 6.0221409e23  # atoms
    extra_ions: int = int(round(n_ions, 0))
    ions: dict[str, int] = {
        "charge_cation_num": simlify_config.solution.charge_cation_extra,
        "charge_anion_num": simlify_config.solution.charge_anion_extra,
    }
    ions: dict[str, int] = {k: v + extra_ions for k, v in ions.items()}
    if simlify_config.solution.charge_neutralize:
        if charge_net < 0:
            ions["charge_cation_num"] += abs(int(charge_net))
        elif charge_net > 0:
            ions["charge_anion_num"] += abs(int(charge_net))
    logger.debug("Ions to add {}", ions)
    return ions
