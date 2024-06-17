import argparse
import os
from abc import ABC, abstractmethod
from collections.abc import Collection

from loguru import logger

from ..utils import get_obj_from_string
from .contexts import SimlifyConfig


class SimPrep(ABC):
    r"""Standardized framework for preparing files for molecular simulation runs."""

    def __init__(self):
        pass

    @staticmethod
    @abstractmethod
    def prepare_sim_config(
        simlify_config: SimlifyConfig,
    ) -> SimlifyConfig:
        r"""Preprocessing and validation of a [simulation context]
        [simulation.prep.SimPrep]. This is unique to each simulation package and can
        be customized based on the context.

        Args:
            simlify_config: Specifies options and parameters.

        Returns:
            Bash commands in a list to run one stage of a simulation.
        """
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def get_stage_input_lines(cls, simlify_config: SimlifyConfig) -> Collection[str]:
        r"""Prepare input file lines for a single stage.

        Args:
            context: Specifies options and parameters.

        Returns:
            Input file lines for a single simulations. The lines do not end in `\n`.
        """
        raise NotImplementedError

    @classmethod
    def prepare_slurm(
        cls,
        file_path: str,
        simlify_config: SimlifyConfig,
        write: bool = True,
    ) -> list[str]:
        r"""Prepare slurm submission script lines.

        Args:
            context: Specifies options and parameters.

        """
        logger.info("Preparing slurm submission script")
        slurm_lines: list[str] = simlify_config.slurm.render()

        slurm_lines = [l + "\n" for l in slurm_lines]
        logger.debug("Slurm script:\n{}", "".join(slurm_lines))
        if write:
            logger.info(f"Writing submissions script at {file_path}")
            with open(file_path, mode="w", encoding="utf-8") as f:
                f.writelines(slurm_lines)
        return slurm_lines

    @classmethod
    @abstractmethod
    def get_stage_run_command(cls, simlify_config: SimlifyConfig) -> Collection[str]:
        r"""Prepare bash command to run a single simulation.

        Args:
            context: Specifies options and parameters.

        Returns:
            Bash commands in a list to run one stage of a simulation.

        **Notes:**

        [`prepare_context`][simulation.amber.prep.AmberSimPrep.prepare_context]
        should be ran before this.
        """
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def prepare_stage(
        cls,
        simlify_config: SimlifyConfig,
        run_commands: list[str] | None = None,
        write: bool = True,
    ) -> tuple[list[str], list[str]]:
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def prepare(cls, simlify_config: SimlifyConfig) -> None:
        """Run all steps to prepare simulations.

        Args:
            simlify_config: Simlify configuration.
        """
        raise NotImplementedError


# pylint: disable-next=too-many-arguments
def run_sim_slurm_prep(
    name_job: str,
    dir_work: str,
    path_run_write: str,
    path_slurm_write: str,
    prep_class_string: str,
    simlify_config: SimlifyConfig,
) -> None:
    r"""Prepare files for simulations using slurm.

    Args:
        name_job: Unique name for this slurm job.
        dir_work: Path to local directory where we will write the simulation files.
        path_run_write: Local path to write a run script with respect to `dir_write`.
        path_slurm_write: Local path to write a slurm submission script with respect to
            `dir_write`.
        prep_class_string: Import string to a simulation preparation class. For example,
            [`"simlify.simulation.amber.prep.AmberSimPrep"`]
            [simulation.amber.prep.AmberSimPrep].
        simlify_config: Simlify configuration.
    """
    simlify_config.rendering.dir_work = dir_work
    simlify_config.path_slurm_write = os.path.join(dir_work, path_slurm_write)
    simlify_config.path_run_write = os.path.join(dir_work, path_run_write)

    sbatch_options = simlify_config.sbatch_options
    if sbatch_options is not None:
        sbatch_options["job-name"] = name_job
        simlify_config.sbatch_options = sbatch_options
    else:
        raise RuntimeError("sbatch options cannot be None when preparing simulations")

    prep_cls = get_obj_from_string(prep_class_string)
    prep_cls.prepare_sbatch_lines(  # type: ignore
        simlify_config, write=simlify_config.write
    )
    prep_cls.prepare(simlify_config)  # type: ignore


def cli_run_sim_slurm_prep():
    r"""Command-line interface to prepare files for simulations using slurm."""
    parser = argparse.ArgumentParser(
        description="Prepare files for simulations using slurm."
    )
    parser.add_argument(
        "name_job",
        type=str,
        nargs="?",
        help="Unique name for this slurm job.",
    )
    parser.add_argument(
        "dir_write",
        type=str,
        nargs="?",
        help="Directory to write input files.",
    )
    parser.add_argument(
        "path_run_write",
        type=str,
        nargs="?",
        help="Local path to write a run script.",
    )
    parser.add_argument(
        "path_slurm_write",
        type=str,
        nargs="?",
        help="Local path to write a slurm submission script.",
    )
    parser.add_argument(
        "prep_class_string",
        type=str,
        nargs="?",
        help="Import string to a simulation preparation class.",
    )
    parser.add_argument(
        "--yaml",
        type=str,
        nargs="*",
        help="Paths to YAML files to use in decreasing precedence.",
    )
    args = parser.parse_args()

    simlify_config = SimlifyConfig()
    if args.yaml is None:
        args.yaml = []
    for yaml_path in reversed(args.yaml):
        simlify_config.from_yaml(yaml_path)
    run_sim_slurm_prep(
        args.name_job,
        args.dir_write,
        args.path_run_write,
        args.path_slurm_write,
        args.prep_class_string,
        simlify_config,
    )
