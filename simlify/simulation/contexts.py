from typing import Any

import argparse
from collections.abc import Collection, Iterable

from loguru import logger
from ruamel.yaml import YAML

from ..utils import get_obj_from_string


# pylint: disable-next=too-many-instance-attributes
class SimContextManager:
    """Contexts for setting up molecular simulations."""

    # pylint: disable-next=too-many-statements
    def __init__(
        self, yaml_paths: str | Iterable[str] | None = None, **kwargs: dict[str, Any]
    ) -> None:
        """
        Args:
            yaml_paths: Path(s) to YAML file(s) to load into the context.
        """
        # Default values in alphabetical order.
        self.charge_anion_extra: int = 0
        """Number of extra anions of type [`charge_anion_identity`]
        [simulation.contexts.SimContextManager.charge_anion_identity] to add to
        the system. This does not include any ions added if [`charge_neutralize`]
        [simulation.contexts.SimContextManager.charge_neutralize] is `True`.
        """
        self.charge_anion_identity: str = "Cl-"
        """Many simulations include anions to either neutralize charge or prepare the
        solvent environment to have a specific ionic concentration. This specifies the
        anion to add to accomplish this.
        """
        self.charge_cation_extra: int = 0
        """Number of extra cations of type [`charge_cation_identity`]
        [simulation.contexts.SimContextManager.charge_cation_identity] to add to the
        system. This does not include any ions added if [`charge_neutralize`]
        [simulation.contexts.SimContextManager.charge_neutralize] is `True`.
        """
        self.charge_cation_identity: str = "Na+"
        """Many simulations include cations to either neutralize charge or prepare the
        solvent environment to have a specific ionic concentration. This specifies the
        cation to add to accomplish this.
        """
        self.charge_net: int = 0
        """Net charge of the molecular system before any preparation."""
        self.charge_neutralize: bool = True
        """Flag to determine if system charge should be neutralized by placing
        additional ions of type [`charge_cation_identity`]
        [simulation.contexts.SimContextManager.charge_cation_identity]
        or [`charge_anion_identity`]
        [simulation.contexts.SimContextManager.charge_anion_identity] based on
        the value of [`charge_net`]
        [simulation.contexts.SimContextManager.charge_net]
        .
        """
        self.compute_platform: str = "mpi"
        """Desired computational platform to run simulations on.

        **Options:**

        -   `mpi`: Message passing interface for central processing units (CPUs).
        -   `cuda`: Compute Unified Device Architecture (CUDA) for graphics processing
            units (GPUs).
        """
        self.cpu_cores: int | None = None
        """Number of CPU cores to use if [compute_platform]
        [simulation.contexts.SimContextManager.compute_platform] is `mpi`."""
        self.dir_input: str | None = None
        """Path to the directory that contains input files when running the simulation.
        """
        self.dir_input_write: str = ""
        """Directory to write input files with respect to `dir_write` during
        preparation.
        """
        self.dir_output: str | None = None
        """Path to the directory that the simulation will store output files."""
        self.dir_scratch: str | None = None
        """Specify path for scratch directory if desired. If `None`, we do not use
        scratch."""
        self.dir_work: str | None = None
        """Directory to be in when running the simulation."""
        self.dir_write: str = "."
        """Local directory to write input files when preparing simulations."""
        self.extra_lines_topo_gen: Iterable[str] | None = None
        """Extra lines to include when generating a topology."""
        self.ff_dna: str | None = None
        """Molecular mechanics force fields for DNA."""
        self.ff_glycam: str | None = None
        """Molecular mechanics force fields for sugars"""
        self.ff_ions: str | None = None
        """Molecular mechanics force fields for ions."""
        self.ff_lipid: str | None = None
        """Molecular mechanics force fields for lipids."""
        self.ff_protein: str | None = None
        """Molecular mechanics force field used to describe polypeptides."""
        self.ff_rna: str | None = None
        """Molecular mechanics force fields for RNA."""
        self.ff_small_molecule: str | None = None
        """Molecular mechanics force fields for small molecules."""
        self.ff_water: str | None = None
        """Molecular mechanics force fields for water."""
        self.input_kwargs: dict[str, Any] | None = None
        """Simulation keyword arguments for input files."""
        self.name_stage: str | None = None
        """Name or label for simulation stage."""
        self.path_coord: str | None = None
        """File path that this simulation will write coordinates data to.

        ### Amber

        An example path would be `$SLURM_SUBMIT_DIR/outputs/05-prod.nc`; Amber would
        then write the coordinates of the simulation to this file.
        """
        self.path_coord_write: str | None = None
        """Path to write coordinate file during preparation."""
        self.path_coord_prev: str | None = None
        """Path to coordinate file of previous stage."""
        self.path_coord_ref: str | None = None
        """Path to reference coordinate file. This is often used for enforcing
        restraints.
        """
        self.path_input: str | None = None
        """Path to the input file during the simulation."""
        self.path_output: str | None = None
        """Path to the output file during the simulation."""
        self.path_restart: str | None = None
        """Path to write the restart file during the simulation."""
        self.path_restart_prev: str | None = None
        """Path to restart file from previous stage or coordinates to start a
        simulation from.
        """
        self.path_run_write: str | None = None
        """Path to write a run script when preparing simulations."""
        self.path_slurm_write: str | None = None
        """Path to write a slurm submission script when preparing simulations."""
        self.path_topo: str | None = None
        """Path to the topology file when the simulation is running."""
        self.path_topo_write: str | None = None
        """Path to write the topology file during preparation."""
        self.sbatch_options: dict[str, Any] | None = None
        """[`sbatch` options](https://slurm.schedmd.com/sbatch.html#SECTION_OPTIONS)
        for a [slurm](https://slurm.schedmd.com/) submission script.
        Some common options are:
        [job-name](https://slurm.schedmd.com/sbatch.html#OPT_job-name),
        [nodes](https://slurm.schedmd.com/sbatch.html#OPT_nodes),
        [ntasks-per-node](https://slurm.schedmd.com/sbatch.html#OPT_ntasks-per-node),
        [cpus-per-task](https://slurm.schedmd.com/sbatch.html#OPT_cpus-per-task),
        [gpus](https://slurm.schedmd.com/sbatch.html#OPT_gpus),
        [gres](https://slurm.schedmd.com/sbatch.html#OPT_gres),
        [cpus-per-gpu](https://slurm.schedmd.com/sbatch.html#OPT_cpus-per-gpu),
        [chdir](https://slurm.schedmd.com/sbatch.html#OPT_chdir),
        [output](https://slurm.schedmd.com/sbatch.html#OPT_output),
        [error](https://slurm.schedmd.com/sbatch.html#OPT_error),
        [time](https://slurm.schedmd.com/sbatch.html#OPT_time),
        [clusters](https://slurm.schedmd.com/sbatch.html#OPT_clusters),
        [partition](https://slurm.schedmd.com/sbatch.html#OPT_partition),
        [account](https://slurm.schedmd.com/sbatch.html#OPT_account).

        These options are written in the format of `#SBATCH --{key}={value}`.
        """
        self.sbatch_lines: list[str] | None = None
        """Lines for a slurm submission script after [sbatch_options]
        [simulation.contexts.SimContextManager.sbatch_options]
        """
        self.solvent_ionic_strength: float = 0.150
        """Ionic strength of the solvent in mole/L."""
        self.solvent_padding: float = 10.0
        """Padding between solute and box edge to fill with solvent in Angstroms."""
        self.splits: int = 1
        """Split simulation stage into several chunks."""
        self.stages: Collection[dict[str, Any]] | None = None
        """Contexts for successive stages. Stage $i > 0$ is assumed to be restarted from
        stage $i - 1$."""
        self.submit: bool = False
        """Submit and run the simulation."""
        self.verbosity: int | str | None = None
        """Verbosity level for logging."""
        self.write: bool = True
        """Flag to write the files during preparation instead of a dry run."""

        if isinstance(yaml_paths, str):
            yaml_paths = [yaml_paths]
        if yaml_paths is not None:
            for yaml_path in yaml_paths:
                self.from_yaml(yaml_path)
        self.update(kwargs)

    def from_yaml(self, yaml_path: str | None) -> None:
        """Load context information from a YAML file. This will only update data
        contained in the YAML file.

        Args:
            yaml_path: Path to YAML file to load.
        """
        if yaml_path is not None:
            logger.info("Loading YAML context from {}", yaml_path)
            yaml = YAML(typ="safe")
            with open(yaml_path, "r", encoding="utf-8") as f:
                yaml_data = yaml.load(f)
            logger.debug("YAML data:\n{}", yaml_data)
            self.update(yaml_data)
        self.yaml_path = yaml_path

    def update(self, attr_dict: dict[str, Any]) -> None:
        """Update attributes with values from the provided dictionary.

        Args:
            attr_dict: Dictionary containing attribute names and their
            corresponding values.
        """
        logger.debug("Updating context:\n{}", attr_dict)
        for key, value in attr_dict.items():
            setattr(self, key, value)

    def get(self) -> dict[str, Any]:
        """Retrieve the context.

        Returns:
            A dictionary representing the current context.
        """
        # The following line filters methods and attributes like __dict__.
        context = {
            k: v for k, v in vars(self).items() if not callable(v) and "__" not in k
        }
        logger.debug("Retrieved context:\n{}", context)
        return context

    def __enter__(self) -> dict[str, Any]:
        """Enter the context and return the current context as a dictionary."""
        return self.get()

    def __exit__(self, exc_type, exc_value, exc_tb):
        """Exit the context.

        Args:
            exc_type: Type of the exception.
            exc_value: Value of the exception.
            exc_tb: Traceback information.
        """


# pylint: disable-next=too-few-public-methods
class ContextValidator:
    """Base class for validating simulation contexts."""

    # pylint: disable=unused-argument

    @classmethod
    def validate(cls, context_manager: SimContextManager) -> bool:
        """Validate contexts for simulations.

        Args:
            context_manager: A simulation context manager to validate.

        Returns:
            If the context is valid.
        """
        logger.info("Validating with: {}", cls.__name__)
        is_valid = True
        context = context_manager.get()
        for key, value in context.items():
            try:
                checker = getattr(cls, key)
            except AttributeError:
                logger.debug("Cannot check {}", key)
                continue
            logger.debug("Checking {}", key)
            if value is not None:
                if isinstance(checker, tuple):
                    if value not in checker:
                        logger.error("  Invalid: {}", value)
                        is_valid = False
                    else:
                        logger.debug("  Valid: {}", value)
                if callable(checker):
                    is_value_valid = checker(value, context)
                    if not is_value_valid:
                        logger.error("  Invalid: {}", value)
                        is_valid = False
                    else:
                        logger.debug("  Valid: {}", value)
            else:
                logger.debug("  Skipping: None")
        return is_valid

    @staticmethod
    def write(value: Any, context: dict[str, Any]) -> bool:
        """Validate `write`"""
        return True

    @staticmethod
    def extra_cations(value: Any, context: dict[str, Any]) -> bool:
        """Validate `extra_cations`"""
        if not isinstance(value, int):
            logger.error("extra_cations must be `int` type")
            return False
        if value < 0:
            logger.error("extra_cations cannot be negative")
            return False
        return True

    @staticmethod
    def extra_anions(value: Any, context: dict[str, Any]) -> bool:
        """Validate `extra_anions`"""
        if not isinstance(value, int):
            logger.error("extra_anions must be `int` type")
            return False
        if value < 0:
            logger.error("extra_anions cannot be negative")
            return False
        return True


def run_context_yaml_validator(
    yaml_paths: str | Iterable[str], validator_obj_string: str
) -> bool:
    """Validate YAML context.

    Args:
        yaml_paths: Path to YAML file.
        validator_obj_string: String to validator object.

    Returns:
        If the YAML context is valid.
    """
    logger.info("Validating context built from {}", yaml_paths)
    validator_cls = get_obj_from_string(validator_obj_string)
    context_manager = SimContextManager(yaml_paths=yaml_paths)
    is_valid: bool = validator_cls.validate(context_manager)  # type: ignore
    valid_string = "IS"
    if not is_valid:
        valid_string += " NOT"
    logger.info("Context {} valid", valid_string)
    return is_valid


def cli_validate_yaml_context() -> None:
    """Command-line interface for validating YAML context files."""
    parser = argparse.ArgumentParser(description="Validate YAML context file.")
    parser.add_argument(
        "validator_string",
        type=str,
        nargs="?",
        help="Import string for validator class.",
    )
    parser.add_argument(
        "--yaml",
        type=str,
        nargs="+",
        help="Path to YAML files to build a context and validate.",
    )
    args = parser.parse_args()
    is_valid = run_context_yaml_validator(args.yaml, args.validator_string)
    if not is_valid:
        raise RuntimeError("Context is not valid")
