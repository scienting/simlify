import os

from atomea.schemas.workflow.amber.schema import AmberSchemaBase
from loguru import logger

from ..contexts import SimlifyConfig
from ..prep import SimPrep


class AmberSimPrep(SimPrep):
    r"""Prepares files for Amber simulations"""

    @staticmethod
    def prepare_sim_config(
        simlify_config: SimlifyConfig,
    ) -> SimlifyConfig:
        r"""Preprocessing and validating of context for Amber simulations.

        Args:
            simlify_config: Simlify configuration.

        """
        logger.info("Preparing context for Amber simulation")

        if not issubclass(simlify_config.engine, AmberSchemaBase):
            raise TypeError(
                "SimlifyConfig.engine should be a subclass of AmberSchemaBase."
            )

        # Creates all amber-specific paths we will need.
        if simlify_config.label is None:
            raise ValueError("simlify_config.label cannot be None")
        simlify_config.engine.cli.mdin = os.path.join(
            simlify_config.runtime.dir_input, simlify_config.label + ".in"
        )
        simlify_config.engine.cli.mdout = os.path.join(
            simlify_config.runtime.dir_run, simlify_config.label + ".out"
        )
        simlify_config.engine.cli.restrt = os.path.join(
            simlify_config.runtime.dir_run, simlify_config.label + ".rst"
        )
        simlify_config.engine.cli.inpcrd = os.path.join(
            simlify_config.runtime.dir_run, simlify_config.label + ".nc"
        )
        simlify_config.engine.cli.mdinfo = os.path.join(
            simlify_config.runtime.dir_run, simlify_config.label + ".mdinfo"
        )

        # Checks to see if we will be writing files to a directory later.
        if simlify_config.dir_write is not None:
            if not os.path.exists(simlify_config.dir_write):
                os.makedirs(simlify_config.dir_write, exist_ok=True)

        if simlify_config.cli.compute_platform == "mpi":
            simlify_config.temp.cpu_cores = int(
                simlify_config.slurm.nodes * simlify_config.slurm.ntasks_per_node
            )

        return simlify_config

    @classmethod
    def get_stage_run_command(cls, simlify_config: SimlifyConfig) -> list[str]:
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
        [`simlify_config`][simulation.contexts.SimlifyConfig].

        -   `name_stage`: Unique label for this simulation stage. This will be used
            to build file paths.
        -   `path_input`: Path to input file for this Amber simulations.
        -   `dir_output`: Path to final output directory.
        -   `dir_scratch`: Path to scratch directory if desired. Will run the
            calculations here and then copy to `dir_output`.
        -   `compute_platform`: Computational platform to run the simulation on.
        -   `cpu_cores`: Number of cores to use for `mpi` simulations if requested.
        """
        use_scratch = bool(simlify_config.dir_scratch is not None)

        stage_commands = ["", f"echo 'Starting {simlify_config.label}'", "date"]

        if simlify_config.runtime.use_scratch:
            # Adds commands to check if split was already ran.
            check_path = os.path.join(
                simlify_config.dir_output, simlify_config.label_stage + ".rst"
            )
            stage_commands = ["    " + line for line in stage_commands]
            stage_commands.insert(0, f"FILE={check_path}")
            stage_commands.insert(1, 'if [ ! -f "$FILE" ]; then')

        amber_command = simlify_config.engine.cli.render()
        if "mpi" in simlify_config.engine.cli.compute_platform.lower():
            amber_command = (
                f"mpirun -np {simlify_config.temp.cpu_cores} " + amber_command
            )

        # TODO: Where do we update with the previous restart file path?
        logger.debug("Amber command: %s", amber_command[:-2])
        stage_commands.append(amber_command)

        if use_scratch:
            # Add indentation to amber command
            stage_commands[-1] = "    " + stage_commands[-1]

            # Adds commands to move scratch files to dir_output.
            stage_commands.extend(
                [
                    f"    mv {simlify_config.engine.cli.mdout} {simlify_config.runtime.dir_output}",
                    f"    mv {simlify_config.engine.cli.restrt} {simlify_config.runtime.dir_output}",
                    f"    mv {simlify_config.engine.cli.inpcrd} {simlify_config.runtime.dir_output}",
                    "fi",
                ]
            )

        return stage_commands

    @classmethod
    def get_stage_input_lines(cls, simlify_config: SimlifyConfig) -> list[str]:
        r"""Prepare input file lines for a single stage.

        Args:
            context: Specifies options and parameters.

        Returns:
            Input file lines for a single simulations. The lines do not end in `\n`.
        """
        logger.info("Preparing input lines for stage {}", simlify_config.label_stage)
        input_lines: list[str] = simlify_config.engine.inputs.render()
        return input_lines

    @classmethod
    def prepare_stage(
        cls,
        simlify_config: SimlifyConfig,
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

        Notes:
            [`prepare_context`][simulation.amber.prep.AmberSimPrep.prepare_context]
            should be ran before this.
        """
        if run_commands is None or len(run_commands) == 0:
            run_commands = ["#!/usr/bin/env bash"]

        # We do not want to change source context in prepare_context, so we do this
        # here.
        if simlify_config.runtime.splits > 1:
            simlify_config.engine.inputs.nstlim = int(
                simlify_config.engine.inputs.nstlim / simlify_config.runtime.splits
            )

        name_stage = simlify_config.label_stage
        for i_split in range(1, simlify_config.runtime.splits + 1):
            if simlify_config.splits > 1:
                name_stage_suffix = f"_split_{i_split:03d}"  # Why n_splits < 1000
            else:
                name_stage_suffix = ""

            simlify_config.label_stage = name_stage + name_stage_suffix
            stage_input_lines = cls.get_stage_input_lines(context)
            if write:
                stage_path_input = os.path.join(
                    simlify_config.dir_write,
                    simlify_config.dir_input_write,
                    simlify_config.label_stage + ".in",
                )
                logger.info("Writing input file at {}", stage_path_input)
                with open(stage_path_input, mode="w", encoding="utf-8") as f:
                    f.writelines([i + "\n" for i in stage_input_lines])

            # Prepare script to run calculations.
            simlify_config.path_coord = os.path.join(
                simlify_config.dir_run, simlify_config.label_stage + ".nc"
            )

            logger.info("Adding stage's amber command to run script")
            stage_commands = cls.get_stage_run_command(context)
            run_commands.extend(stage_commands)
        return stage_input_lines, run_commands

    @classmethod
    def prepare(cls, simlify_config: SimlifyConfig) -> None:
        """Run all steps to prepare simulations.

        Args:
            simlify_config: Simlify configuration.
        """
        logger.info("Prepare simulations")
        multiple_stages: bool = bool(simlify_config.stages is not None)
        if multiple_stages:
            logger.debug("There are multiple stages")
            n_stages: int = len(simlify_config.stages)  # type: ignore
            simlify_config.update(simlify_config.stages[0])  # type: ignore
        else:
            logger.debug("There is one stage.")
            n_stages = 1

        run_commands: list[str] = []
        for i_stage in range(1, n_stages + 1):
            logger.info("Preparing stage {}", i_stage - 1)
            simlify_config = cls.prepare_context(simlify_config)
            context = simlify_config.get()
            _, run_commands = cls.prepare_stage(
                context, run_commands, simlify_config.write
            )

            simlify_config.path_restart_prev = simlify_config.path_restart
            simlify_config.path_coord_prev = simlify_config.path_coord

            if multiple_stages and i_stage < n_stages:
                simlify_config.update(simlify_config.stages[i_stage])  # type: ignore

        if simlify_config.write:
            logger.debug("Writing run script at {}", simlify_config.path_run_write)
            with open(simlify_config.path_run_write, "w", encoding="utf-8") as f:
                f.writelines([l + "\n" for l in run_commands])
