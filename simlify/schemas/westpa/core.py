from pydantic import BaseModel, Field

from simlify.schemas.utils import YamlIO
from simlify.schemas.westpa import WestpaEnv
from simlify.schemas.westpa.cfg.core import WestpaConfigConfig


class WestpaConfig(BaseModel, YamlIO):
    """The root configuration for `west.cfg` files."""

    west: WestpaConfigConfig = Field(default_factory=WestpaConfigConfig)
    env: WestpaEnv = Field(default_factory=WestpaEnv)
