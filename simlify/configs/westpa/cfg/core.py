from pydantic import BaseModel, Field

from simlify.configs.utils import Render
from simlify.configs.westpa.cfg import (
    DataConfig,
    ExecutableConfig,
    PropagationConfig,
    SystemConfig,
)


class WestpaConfigConfig(BaseModel, Render):
    """
    Configuration class for wrapping WESTPA configuration options. These options
    are scattered across a variety of sources shown below. We do our best to aggregate
    and explain each parameter.

    -   [GitHub wiki](https://github.com/westpa/westpa/wiki/Configuration-File),
    -   [readthedocs](https://westpa.readthedocs.io/en/latest/users_guide/west/setup.html#configuration-file),
    """

    system: SystemConfig = Field(default_factory=SystemConfig)
    """
    Configuration options for the physical systems we are simulating.
    """
    propagation: PropagationConfig = Field(default_factory=PropagationConfig)
    """
    How to propagate the simulations.
    """
    data: DataConfig = Field(default_factory=DataConfig)
    """
    Specify what, and how, data should be stored.
    """
    executable: ExecutableConfig = Field(default_factory=ExecutableConfig)
    """
    Provides information for the executable.
    """
