from typing import Any

import os

from loguru import logger

from ..contexts import SimContextManager
from ..prep import SimPrep
from .contexts import AmberContextValidator


class AmberSimPrep(SimPrep):
    r"""Prepares files for Amber simulations"""

    @staticmethod
    def prepare_context(
        sim_context_manager: SimContextManager,
    ) -> SimContextManager:
        r"""Preprocessing and validating of context for Amber simulations.

        Args:
            sim_context_manager: Context manager for simulations.

        """
        logger.info("Preparing context for Amber simulation")
        context = sim_context_manager.get()

        dir_run = context["dir_output"]

        use_scratch = bool(context["dir_scratch"] is not None)
        if use_scratch:
            dir_run = context["dir_scratch"]
        context["dir_run"] = dir_run

        # Creates all amber-specific paths we will need.
        context["path_input"] = os.path.join(
            context["dir_input"], context["name_stage"] + ".in"
        )
        context["path_output"] = os.path.join(dir_run, context["name_stage"] + ".out")
        context["path_restart"] = os.path.join(dir_run, context["name_stage"] + ".rst")
        context["path_coord"] = os.path.join(dir_run, context["name_stage"] + ".nc")
        context["path_mdinfo"] = os.path.join(
            dir_run, context["name_stage"] + ".mdinfo"
        )

        dir_input_write = os.path.join(context["dir_write"], context["dir_input_write"])
        if not os.path.exists(dir_input_write):
            os.makedirs(dir_input_write, exist_ok=True)

        if use_scratch:
            # We have to use absolute paths with scratch to ensure nothing gets
            # overwritten.
            for k, v in context.items():
                if k[-5:] in ("_path", "_dir"):
                    if v is not None:
                        if "$slurm" not in v.lower():
                            context[k] = os.path.abspath(v)

        if "sbatch_options" in context.keys():
            if context["compute_platform"] == "mpi":
                if context["cpu_cores"] is None:
                    n_nodes = context["sbatch_options"]["nodes"]
                    ntasks_per_node = context["sbatch_options"]["ntasks-per-node"]
                    context["cpu_cores"] = int(n_nodes * ntasks_per_node)

        sim_context_manager.update(context)

        is_valid = AmberContextValidator.validate(sim_context_manager)
        if not is_valid:
            raise RuntimeError("Context is not valid for Amber!")
        return sim_context_manager

    @classmethod
    def get_stage_run_command(cls, context: dict[str, Any]) -> list[str]:
        r"""Prepare bash command to run a single Amber simulation using pmemd.

        Args:
            context: Specifies options and parameters.

        Returns:
            Bash commands in a list to run one stage of a simulation.

        **Notes:**

        [`prepare_context`][simulation.amber.prep.AmberSimPrep.prepare_context]
        should be ran before this.

        **Uses:**

        The following attributes are possibly used here and should be specified in
        [`sim_context_manager`][simulation.contexts.SimContextManager].

        -   `name_stage`: Unique label for this simulation stage. This will be used
            to build file paths.
        -   `path_input`: Path to input file for this Amber simulations.
        -   `dir_output`: Path to final output directory.
        -   `dir_scratch`: Path to scratch directory if desired. Will run the
            calculations here and then copy to `dir_output`.
        -   `compute_platform`: Computational platform to run the simulation on.
        -   `cpu_cores`: Number of cores to use for `mpi` simulations if requested.
        """
        use_scratch = bool(context["dir_scratch"] is not None)

        stage_commands = ["", f"echo 'Starting {context['name_stage']}'", "date"]

        if use_scratch:
            # Adds commands to check if split was already ran.
            check_path = os.path.join(
                context["dir_output"], context["name_stage"] + ".rst"
            )
            stage_commands = ["    " + line for line in stage_commands]
            stage_commands.insert(0, f"FILE={check_path}")
            stage_commands.insert(1, 'if [ ! -f "$FILE" ]; then')

        if context["compute_platform"] == "mpi":
            amber_command = f"mpirun -np {context['cpu_cores']} pmemd.MPI "
        elif context["compute_platform"] == "cuda":
            amber_command = "pmemd.cuda "
        amber_command += f"-O -i {context['path_input']} "
        amber_command += (
            f"-o {context['path_output']} -c {context['path_restart_prev']} "
        )
        amber_command += f"-p {context['path_topo']} "
        amber_command += f"-r {context['path_restart']} -x {context['path_coord']} "
        amber_command += (
            f"-ref {context['path_coord_ref']} -inf {context['path_mdinfo']}"
        )
        logger.debug("Amber command: %s", amber_command[:-2])
        stage_commands.append(amber_command)

        if use_scratch:
            # Add indentation to amber command
            stage_commands[-1] = "    " + stage_commands[-1]

            # Adds commands to move scratch files to dir_output.
            stage_commands.extend(
                [
                    f"    mv {context['path_output']} {context['dir_output']}",
                    f"    mv {context['path_restart']} {context['dir_output']}",
                    f"    mv {context['path_coord']} {context['dir_output']}",
                    "fi",
                ]
            )

        return stage_commands

    @classmethod
    def get_stage_input_lines(cls, context: dict[str, Any]) -> list[str]:
        r"""Prepare input file lines for a single stage.

        Args:
            context: Specifies options and parameters.

        Returns:
            Input file lines for a single simulations. The lines do not end in `\n`.
        """
        logger.info("Preparing input lines for stage {}", context["name_stage"])
        input_lines = [context["name_stage"], "&cntrl"]
        for key, value in context["input_kwargs"].items():
            if key in ("restraintmask", "timask1", "scmask1", "timask2", "scmask2"):
                value = f'"{value}"'
            line_to_add = f"    {key}={value},"
            logger.debug("Adding input line: {}", line_to_add.strip())
            input_lines.append(line_to_add)
        input_lines.append("&end")
        return input_lines

    @classmethod
    def prepare_stage(
        cls,
        context: dict[str, Any],
        run_commands: list[str] | None = None,
        write: bool = True,
    ) -> tuple[list[str], list[str]]:
        r"""Write input files for a simulation stage and builds bash commands to
        run all splits.

        Args:
            context: Specifies options and parameters.
            run_commands: Cumulative run commands for all desired stages.
            write: Write input file to disk.

        Returns:
            Input file lines for this stage.

            Updated `run_commands` including this stage.

        **Notes:**

        [`prepare_context`][simulation.amber.prep.AmberSimPrep.prepare_context]
        should be ran before this.
        """
        if run_commands is None or len(run_commands) == 0:
            run_commands = ["#!/usr/bin/env bash"]

        # We do not want to change source context in prepare_context, so we do this
        # here.
        if context["splits"] > 1:
            context["input_kwargs"]["nstlim"] = int(
                context["input_kwargs"]["nstlim"] / context["splits"]
            )

        name_stage = context["name_stage"]
        for i_split in range(1, context["splits"] + 1):
            if context["splits"] > 1:
                name_stage_suffix = f"_split_{i_split:03d}"  # Why n_splits < 1000
            else:
                name_stage_suffix = ""

            context["name_stage"] = name_stage + name_stage_suffix
            stage_input_lines = cls.get_stage_input_lines(context)
            if write:
                stage_path_input = os.path.join(
                    context["dir_write"],
                    context["dir_input_write"],
                    context["name_stage"] + ".in",
                )
                logger.info("Writing input file at {}", stage_path_input)
                with open(stage_path_input, mode="w", encoding="utf-8") as f:
                    f.writelines([i + "\n" for i in stage_input_lines])

            # Prepare script to run calculations.
            context["path_coord"] = os.path.join(
                context["dir_run"], context["name_stage"] + ".nc"
            )

            logger.info("Adding stage's amber command to run script")
            stage_commands = cls.get_stage_run_command(context)
            run_commands.extend(stage_commands)
        return stage_input_lines, run_commands

    @classmethod
    def prepare(cls, sim_context_manager: SimContextManager) -> None:
        """Run all steps to prepare simulations.

        Args:
            sim_context_manager: Context manager for simulations.
        """
        logger.info("Prepare simulations")
        multiple_stages: bool = bool(sim_context_manager.stages is not None)
        if multiple_stages:
            logger.debug("There are multiple stages")
            n_stages: int = len(sim_context_manager.stages)  # type: ignore
            sim_context_manager.update(sim_context_manager.stages[0])  # type: ignore
        else:
            logger.debug("There is one stage.")
            n_stages = 1

        run_commands: list[str] = []
        for i_stage in range(1, n_stages + 1):
            logger.info("Preparing stage {}", i_stage - 1)
            sim_context_manager = cls.prepare_context(sim_context_manager)
            context = sim_context_manager.get()
            _, run_commands = cls.prepare_stage(context, run_commands, context["write"])

            sim_context_manager.path_restart_prev = sim_context_manager.path_restart
            sim_context_manager.path_coord_prev = sim_context_manager.path_coord

            if multiple_stages and i_stage < n_stages:
                sim_context_manager.update(sim_context_manager.stages[i_stage])  # type: ignore

        if context["write"]:
            logger.debug("Writing run script at {}", context["path_run_write"])
            with open(context["path_run_write"], "w", encoding="utf-8") as f:
                f.writelines([l + "\n" for l in run_commands])
