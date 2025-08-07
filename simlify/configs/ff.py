# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

from atomea.schemas import YamlIO
from pydantic import BaseModel

from simlify.configs.utils import YamlIO


class ForcefieldConfigBase(BaseModel, YamlIO):
    dna: str | None = None
    """Molecular mechanics force fields for DNA."""

    glycam: str | None = None
    """Molecular mechanics force fields for sugars"""

    ions: str | None = None
    """Molecular mechanics force fields for ions."""

    lipid: str | None = None
    """Molecular mechanics force fields for lipids."""

    protein: str | None = None
    """Molecular mechanics force field used to describe polypeptides."""

    rna: str | None = None
    """Molecular mechanics force fields for RNA."""

    small_molecule: str | None = None
    """Molecular mechanics force fields for small molecules."""

    water: str | None = None
    """Molecular mechanics force fields for water."""
