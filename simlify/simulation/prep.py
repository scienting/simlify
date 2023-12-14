from typing import Any

import argparse
import os
from abc import ABC, abstractmethod
from collections.abc import Collection

from loguru import logger

from ..utils import get_obj_from_string
from .contexts import SimContextManager


class SimPrep(ABC):
    r"""Standardized framework for preparing files for molecular simulation runs."""

    def __init__(self):
        pass

    @staticmethod
    @abstractmethod
    def prepare_context(
        sim_context_manager: SimContextManager,
    ) -> SimContextManager:
        r"""Preprocessing and validation of a [simulation context]
        [simulation.prep.SimPrep]. This is unique to each simulation package and can
        be customized based on the context.

        Args:
            sim_context_manager: Specifies options and parameters.

        Returns:
            Bash commands in a list to run one stage of a simulation.
        """
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def get_stage_input_lines(cls, context: dict[str, Any]) -> Collection[str]:
        r"""Prepare input file lines for a single stage.

        Args:
            context: Specifies options and parameters.

        Returns:
            Input file lines for a single simulations. The lines do not end in `\n`.
        """
        raise NotImplementedError

    @classmethod
    def prepare_sbatch_lines(
        cls,
        context: dict[str, Any],
        write: bool = True,
    ) -> list[str]:
        r"""Prepare slurm submission script lines.

        Args:
            context: Specifies options and parameters.

        """
        logger.info("Preparing slurm submission script")
        sbatch_lines = ["#!/bin/bash"]
        for k, v in context["sbatch_options"].items():
            sbatch_lines.append(f"#SBATCH --{k}={v}")
        sbatch_lines.extend(context["sbatch_lines"])

        sbatch_lines = [l + "\n" for l in sbatch_lines if isinstance(l, str)]
        logger.debug("Slurm script:\n{}", "".join(sbatch_lines))
        if write:
            logger.info("Writing submissions script at {}", context["path_slurm_write"])
            with open(context["path_slurm_write"], mode="w", encoding="utf-8") as f:
                f.writelines(sbatch_lines)
        return sbatch_lines

    @classmethod
    @abstractmethod
    def get_stage_run_command(cls, context: dict[str, Any]) -> Collection[str]:
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
        context: dict[str, Any],
        run_commands: list[str] | None = None,
        write: bool = True,
    ) -> tuple[list[str], list[str]]:
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def prepare(cls, sim_context_manager: SimContextManager) -> None:
        """Run all steps to prepare simulations.

        Args:
            sim_context_manager: Context manager for simulations.
        """
        raise NotImplementedError


# pylint: disable-next=too-many-arguments
def run_sim_slurm_prep(
    name_job: str,
    dir_write: str,
    path_run_write: str,
    path_slurm_write: str,
    prep_class_string: str,
    sim_context_manager: SimContextManager,
) -> None:
    r"""Prepare files for simulations using slurm.

    Args:
        name_job: Unique name for this slurm job..
        dir_write: Path to local directory where we will write the simulation files.
        path_run_write: Local path to write a run script with respect to `dir_write`.
        path_slurm_write: Local path to write a slurm submission script with respect to
            `dir_write`.
        prep_class_string: Import string to a simulation preparation class. For example,
            [`"simlify.simulation.amber.prep.AmberSimPrep"`]
            [simulation.amber.prep.AmberSimPrep].
        sim_context_manager: Context manager for simulations.
    """
    sim_context_manager.dir_write = dir_write
    sim_context_manager.path_slurm_write = os.path.join(dir_write, path_slurm_write)
    sim_context_manager.path_run_write = os.path.join(dir_write, path_run_write)

    sbatch_options = sim_context_manager.sbatch_options
    if sbatch_options is not None:
        sbatch_options["job-name"] = name_job
        sim_context_manager.sbatch_options = sbatch_options
    else:
        raise RuntimeError("sbatch options cannot be None when preparing simulations")

    prep_cls = get_obj_from_string(prep_class_string)
    prep_cls.prepare_sbatch_lines(  # type: ignore
        sim_context_manager.get(), write=sim_context_manager.write
    )
    prep_cls.prepare(sim_context_manager)  # type: ignore


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

    sim_context_manager = SimContextManager()
    if args.yaml is None:
        args.yaml = []
    for yaml_path in reversed(args.yaml):
        sim_context_manager.from_yaml(yaml_path)
    run_sim_slurm_prep(
        args.name_job,
        args.dir_write,
        args.path_run_write,
        args.path_slurm_write,
        args.prep_class_string,
        sim_context_manager,
    )
