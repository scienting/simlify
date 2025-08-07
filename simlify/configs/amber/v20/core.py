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
from simlify.configs.amber.v20 import (
    Amber20CLI,
    Amber20Forcefield,
    Amber20Inputs,
)


class Amber20Config(AmberConfigBase):
    r"""Amber 20 config for simulation contexts."""

    inputs: Amber20Inputs = Field(default_factory=Amber20Inputs)

    cli: Amber20CLI = Field(default_factory=Amber20CLI)

    ff: Amber20Forcefield = Field(default_factory=Amber20Forcefield)
