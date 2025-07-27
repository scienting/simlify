# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

from atomea.schemas.io import YamlIO
from pydantic import BaseModel


class SolutionConfig(BaseModel, YamlIO):
    charge_anion_extra: int = 0
    """Number of extra anions of type [`charge_anion_identity`]
    [schemas.SolutionConfig.charge_anion_identity] to add to
    the system. This does not include any ions added if [`charge_neutralize`]
    [schemas.SolutionConfig.charge_neutralize] is `True`.
    """

    charge_anion_identity: str = "Cl-"
    """Many simulations include anions to either neutralize charge or prepare the
    solvent environment to have a specific ionic concentration. This specifies the
    anion to add to accomplish this.
    """

    charge_cation_extra: int = 0
    """Number of extra cations of type [`charge_cation_identity`]
    [schemas.SolutionConfig.charge_cation_identity] to add to the
    system. This does not include any ions added if [`charge_neutralize`]
    [schemas.SolutionConfig.charge_neutralize] is `True`.
    """

    charge_cation_identity: str = "Na+"
    """Many simulations include cations to either neutralize charge or prepare the
    solvent environment to have a specific ionic concentration. This specifies the
    cation to add to accomplish this.
    """

    charge_net: int = 0
    """Net charge of the molecular system before any preparation."""

    charge_neutralize: bool = True
    """Flag to determine if system charge should be neutralized by placing
    additional ions of type [`charge_cation_identity`]
    [schemas.SolutionConfig.charge_cation_identity]
    or [`charge_anion_identity`]
    [schemas.SolutionConfig.charge_anion_identity] based on
    the value of [`charge_net`]
    [schemas.SolutionConfig.charge_net].
    """

    solvent_ionic_strength: float = 0.150
    """Ionic strength of the solvent in mol/L."""

    solvent_padding: float = 10.0
    """Padding between solute and box edge to fill with solvent in Angstroms."""
