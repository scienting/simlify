from collections.abc import Iterable

from atomea.schemas.io import YamlIO
from atomea.schemas.workflow.slurm import SlurmSchema
from pydantic import BaseModel, Field


class ForcefieldConfig(BaseModel, YamlIO):

    dna: str | None = None
    """Molecular mechanics force fields for DNA."""

    glycam: str | None = None
    """Molecular mechanics force fields for sugars"""

    ions: str | None = None
    """Molecular mechanics force fields for ions."""

    lipid: str | None = None
    """Molecular mechanics force fields for lipids."""

    protein: str | None = None
    """Molecular mechanics force field used to describe polypeptides."""

    rna: str | None = None
    """Molecular mechanics force fields for RNA."""

    small_molecule: str | None = None
    """Molecular mechanics force fields for small molecules."""

    water: str | None = None
    """Molecular mechanics force fields for water."""


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


class SimlifyConfig(BaseModel, YamlIO):
    """Contexts for setting up molecular simulations."""

    ff: ForcefieldConfig = Field(default_factory=ForcefieldConfig)

    slurm: SlurmSchema = Field(default_factory=SlurmSchema)

    solution: SolutionConfig = Field(default_factory=SolutionConfig)

    dir_input: str = "."
    """Path to the directory that contains input files when running the simulation.
    """

    dir_output: str = "."
    """Path to the directory that the simulation will store output files."""

    dir_scratch: str | None = None
    """Specify path for scratch directory if desired. If `None`, we do not use
    scratch."""

    dir_work: str = "."
    """Directory to be in when running the simulation."""

    dir_write: str = "."
    """Local directory to write input files when preparing simulations."""

    extra_lines_topo_gen: Iterable[str] | None = None
    """Extra lines to include when generating a topology."""

    name_stage: str | None = None
    """Name or label for simulation stage."""

    splits: int = 1
    """Split simulation stage into several chunks."""
