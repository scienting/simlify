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
from simlify.configs.amber.v22 import (
    Amber22CLI,
    Amber22Forcefield,
    Amber22Inputs,
)


class Amber22Config(AmberConfigBase):
    r"""Amber 22 config for simulation contexts."""

    inputs: Amber22Inputs = Field(default_factory=Amber22Inputs)

    cli: Amber22CLI = Field(default_factory=Amber22CLI)

    ff: Amber22Forcefield = Field(default_factory=Amber22Forcefield)
