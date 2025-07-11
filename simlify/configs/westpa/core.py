from pydantic import BaseModel, Field

from simlify.configs.utils import YamlIO
from simlify.configs.westpa import WestpaEnv
from simlify.configs.westpa.cfg.core import WestpaConfigConfig


class WestpaConfig(BaseModel, YamlIO):
    """The root configuration for `west.cfg` files."""

    west: WestpaConfigConfig = Field(default_factory=WestpaConfigConfig)
    env: WestpaEnv = Field(default_factory=WestpaEnv)
