from typing import Any

from atomea.schemas.io import YamlIo
from pydantic import BaseModel


class GninaSchemaBase(BaseModel, YamlIo):
    """Base model for the Gnina schema."""

    inputs: Any = NotImplementedError

    cli: Any = NotImplementedError
