# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

from pydantic import Field

from simlify.configs.amber import AmberConfigBase
from simlify.configs.amber.v18 import (
    Amber18CLI,
    Amber18Forcefield,
    Amber18Inputs,
)


class Amber18Config(AmberConfigBase):
    r"""Amber 18 config for simulation contexts."""

    inputs: Amber18Inputs = Field(default_factory=Amber18Inputs)

    cli: Amber18CLI = Field(default_factory=Amber18CLI)

    ff: Amber18Forcefield = Field(default_factory=Amber18Forcefield)
