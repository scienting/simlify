# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

from pydantic import Field

from simlify.schemas.amber import AmberSchemaBase
from simlify.schemas.amber.v22 import (
    Amber22CLI,
    Amber22Forcefield,
    Amber22Inputs,
)


class Amber22Schema(AmberSchemaBase):
    r"""Amber 22 schema for simulation contexts."""

    inputs: Amber22Inputs = Field(default_factory=Amber22Inputs)

    cli: Amber22CLI = Field(default_factory=Amber22CLI)

    ff: Amber22Forcefield = Field(default_factory=Amber22Forcefield)
