from typing import Any

from pydantic import BaseModel

from simlify.configs.utils import YamlIO


class AmberConfigBase(BaseModel, YamlIO):
    r"""Validate Amber contexts."""

    inputs: Any = NotImplementedError

    cli: Any = NotImplementedError

    ff: Any = NotImplementedError
