from atomea.schemas.io import YamlIO
from pydantic import BaseModel


class TopologyConfig(BaseModel, YamlIO):
    """Topology configuration."""

    append_lines: list[str] = []
    """Extra lines to include when generating a topology."""
