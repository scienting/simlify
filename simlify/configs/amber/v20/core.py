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
