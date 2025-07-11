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
