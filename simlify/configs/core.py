# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

from typing import Any

from pydantic import BaseModel, Field

from simlify.configs import (
    RunConfig,
    SlurmConfig,
    SolutionConfig,
    TopologyConfig,
)
from simlify.configs.utils import Render, YamlIO


class SimlifyConfig(BaseModel, YamlIO, Render):
    """Contexts for setting up molecular simulations.

    This class aggregates various configuration required to define and
    execute molecular simulations using the Simlify workflow. It includes
    settings for job submission (Slurm), solution conditions, topology generation,
    runtime parameters, and the specific simulation engine to be used.
    """

    slurm: SlurmConfig = Field(default_factory=SlurmConfig)
    """
    Configuration for submitting jobs to a Slurm workload manager.
    This includes parameters such as the number of nodes, tasks per node,
    partition, time limit, and other Slurm-specific settings.
    """

    solution: SolutionConfig = Field(default_factory=SolutionConfig)
    """
    Configuration for defining the solution environment of the
    molecular system. This includes parameters related to solvent, ions,
    and charge neutralization.
    """

    topology: TopologyConfig = Field(default_factory=TopologyConfig)
    """
    Configuration for generating the topology and parameter files
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
    Configuration for specifying parameters related to the simulation
    runtime environment. This includes working directories, input/output paths,
    and options for splitting the simulation into chunks.
    """

    engine: Any = None
    """
    Workflow for the molecular simulation engine (e.g.,
    [`Amber22Config`][simlify.configs.amber.v22.Amber22Config]).
    This attribute holds the specific engine configuration, allowing Simlify
    to interact with different simulation packages.
    """

    label: str | None = None
    """
    A descriptive label for this specific simulation configuration. This can
    be useful for identifying and organizing different simulation setups.
    """
