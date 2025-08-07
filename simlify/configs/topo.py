# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

from atomea.schemas.io import YamlIO
from pydantic import BaseModel

from simlify.configs.utils import YamlIO


class TopologyConfig(BaseModel, YamlIO):
    """Topology configuration."""

    append_lines: list[str] = []
    """Extra lines to include when generating a topology."""
