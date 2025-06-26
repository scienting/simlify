from pydantic import BaseModel

from simlify.schemas.utils import YamlIO


class TopologyConfig(BaseModel, YamlIO):
    """Topology configuration."""

    append_lines: list[str] = []
    """Extra lines to include when generating a topology."""
