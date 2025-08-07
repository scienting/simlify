from typing import Any

from pydantic import BaseModel

from simlify.configs.utils import YamlIO


class GninaConfigBase(BaseModel, YamlIO):
    r"""Gnina config."""
