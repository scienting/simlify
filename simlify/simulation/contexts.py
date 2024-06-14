from typing import Any

from collections.abc import Iterable

from atomea.schemas.io import YamlIO
from atomea.schemas.workflow.slurm import SlurmSchema
from pydantic import BaseModel, Field


class SolutionConfig(BaseModel, YamlIO):
    charge_anion_extra: int = 0
    """Number of extra anions of type [`charge_anion_identity`]
    [simulation.contexts.SimlifyConfig.charge_anion_identity] to add to
    the system. This does not include any ions added if [`charge_neutralize`]
    [simulation.contexts.SimlifyConfig.charge_neutralize] is `True`.
    """

    charge_anion_identity: str = "Cl-"
    """Many simulations include anions to either neutralize charge or prepare the
    solvent environment to have a specific ionic concentration. This specifies the
    anion to add to accomplish this.
    """

    charge_cation_extra: int = 0
    """Number of extra cations of type [`charge_cation_identity`]
    [simulation.contexts.SimlifyConfig.charge_cation_identity] to add to the
    system. This does not include any ions added if [`charge_neutralize`]
    [simulation.contexts.SimlifyConfig.charge_neutralize] is `True`.
    """

    charge_cation_identity: str = "Na+"
    """Many simulations include cations to either neutralize charge or prepare the
    solvent environment to have a specific ionic concentration. This specifies the
    cation to add to accomplish this.
    """

    charge_net: int = 0
    """Net charge of the molecular system before any preparation."""

    charge_neutralize: bool = True
    """Flag to determine if system charge should be neutralized by placing
    additional ions of type [`charge_cation_identity`]
    [simulation.contexts.SimlifyConfig.charge_cation_identity]
    or [`charge_anion_identity`]
    [simulation.contexts.SimlifyConfig.charge_anion_identity] based on
    the value of [`charge_net`]
    [simulation.contexts.SimlifyConfig.charge_net].
    """

    solvent_ionic_strength: float = 0.150
    """Ionic strength of the solvent in mole/L."""

    solvent_padding: float = 10.0
    """Padding between solute and box edge to fill with solvent in Angstroms."""


class TempSchema(BaseModel, YamlIO):
    """Provides a pydantic model to put temporary information."""


class RenderingConfig(BaseModel, YamlIO):
    """Configuration for rendering files to write."""

    dir_work: str = "."
    """
    Working directory to write files.
    """

    dir_input: str = "."
    """
    Path to a directory relative to
    [`dir_work`][simulation.contexts.RenderingConfig.dir_work] that will contain
    input files.
    """

    dir_output: str = "."
    """
    Path to a directory relative to
    [`dir_work`][simulation.contexts.RenderingConfig.dir_work] that the simulation will
    store output files.
    """


class RuntimeConfig(BaseModel, YamlIO):
    """Configuration options during the simulation runtime."""

    dir_work: str = "."
    """
    Working directory during runtime. This can be the current work directory or
    a workload manager directory such as `$SLURM_SUBMIT_DIR`.
    """

    dir_run: str = "."
    """
    Directory that calculations are performed in. This can be the same as
    [`dir_work`][simulation.contexts.RuntimeConfig.dir_work]
    or some scratch space like `$SLURM_SCRATCH`.
    """

    dir_input: str = "."
    """
    Path to a directory relative to
    [`dir_work`][simulation.contexts.RuntimeConfig.dir_work] that will contain
    input files.
    """

    dir_output: str = "."
    """
    Path to a directory relative to
    [`dir_work`][simulation.contexts.RuntimeConfig.dir_work] that the simulation will
    store output files.
    """

    splits: int = 1
    """Split simulation stage into several chunks."""


class SimlifyConfig(BaseModel, YamlIO):
    """Contexts for setting up molecular simulations."""

    slurm: SlurmSchema = Field(default_factory=SlurmSchema)

    solution: SolutionConfig = Field(default_factory=SolutionConfig)

    temp: TempSchema = Field(default_factory=TempSchema, exclude=True)

    rendering: RenderingConfig = Field(default_factory=RenderingConfig)

    runtime: RuntimeConfig = Field(default_factory=RuntimeConfig)

    engine: Any = None
    """
    Atomea workflow schema for the molecular simulation engine (e.g.,
    `Amber22Schema`).
    """

    extra_lines_topo_gen: Iterable[str] | None = None
    """Extra lines to include when generating a topology."""

    label: str | None = None
    """Label for this specific simulation."""
