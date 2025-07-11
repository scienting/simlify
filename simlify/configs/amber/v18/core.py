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
