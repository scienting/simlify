from typing import Any

from pydantic import BaseModel, Field

from simlify.schemas import (
    RunConfig,
    SlurmSchema,
    SolutionConfig,
    TopologyConfig,
)
from simlify.schemas.utils import Render, YamlIO


class SimlifyConfig(BaseModel, YamlIO, Render):
    """Contexts for setting up molecular simulations.

    This class aggregates various configuration schemas required to define and
    execute molecular simulations using the Simlify workflow. It includes
    settings for job submission (Slurm), solution conditions, topology generation,
    runtime parameters, and the specific simulation engine to be used.
    """

    slurm: SlurmSchema = Field(default_factory=SlurmSchema)
    """
    Configuration schema for submitting jobs to a Slurm workload manager.
    This includes parameters such as the number of nodes, tasks per node,
    partition, time limit, and other Slurm-specific settings.
    """

    solution: SolutionConfig = Field(default_factory=SolutionConfig)
    """
    Configuration schema for defining the solution environment of the
    molecular system. This includes parameters related to solvent, ions,
    and charge neutralization.
    """

    topology: TopologyConfig = Field(default_factory=TopologyConfig)
    """
    Configuration schema for generating the topology and parameter files
    required for the molecular simulation. This can include options for
    appending extra lines to the topology file.
    """

    temp: dict[str, Any] = Field(default={}, exclude=True)
    """
    A temporary dictionary for storing intermediate data during the workflow.
    This attribute is excluded from serialization and rendering.
    """

    run: RunConfig = Field(default_factory=RunConfig)
    """
    Configuration schema for specifying parameters related to the simulation
    runtime environment. This includes working directories, input/output paths,
    and options for splitting the simulation into chunks.
    """

    engine: Any = None
    """
    Workflow schema for the molecular simulation engine (e.g.,
    [`Amber22Schema`][simlify.schemas.amber.v22.Amber22Schema]).
    This attribute holds the specific engine configuration, allowing Simlify
    to interact with different simulation packages.
    """

    label: str | None = None
    """
    A descriptive label for this specific simulation configuration. This can
    be useful for identifying and organizing different simulation setups.
    """
