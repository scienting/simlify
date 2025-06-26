from typing import Any

from pydantic import BaseModel

from simlify.schemas.utils import YamlIO


class AmberSchemaBase(BaseModel, YamlIO):
    r"""Validate Amber contexts."""

    inputs: Any = NotImplementedError

    cli: Any = NotImplementedError

    ff: Any = NotImplementedError
