from typing import Any

from atomea.schemas.io import YamlIO
from pydantic import BaseModel


class AmberSchemaBase(BaseModel, YamlIO):
    r"""Validate Amber contexts."""

    inputs: Any = NotImplementedError

    cli: Any = NotImplementedError

    ff: Any = NotImplementedError
