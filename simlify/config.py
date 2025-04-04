from typing import Any

from collections.abc import Iterable

from atomea.schemas import Render, YamlIO
from atomea.schemas.workflow import (
    RunConfig,
    SlurmSchema,
    SolutionConfig,
    TopologyConfig,
)
from pydantic import BaseModel, Field


class SimlifyConfig(BaseModel, YamlIO, Render):
    """Contexts for setting up molecular simulations."""

    slurm: SlurmSchema = Field(default_factory=SlurmSchema)

    solution: SolutionConfig = Field(default_factory=SolutionConfig)

    topology: TopologyConfig = Field(default_factory=TopologyConfig)

    temp: dict[str, Any] = Field(default={}, exclude=True)

    run: RunConfig = Field(default_factory=RunConfig)

    engine: Any = None
    """
    Atomea workflow schema for the molecular simulation engine (e.g.,
    `Amber22Schema`).
    """

    extra_lines_topo_gen: Iterable[str] | None = None
    """Extra lines to include when generating a topology."""

    label: str | None = None
    """Label for this specific simulation."""
