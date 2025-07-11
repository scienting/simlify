import os

from loguru import logger

from simlify import SimlifyConfig
from simlify.simulation.prep import SimPrep


class AmberSimPrep(SimPrep):
    r"""Prepares files for Amber molecular dynamics simulations.

    This class implements the `SimPrep` abstract base class to provide a concrete
    workflow for preparing input files, run scripts, and other necessary files
    specifically for running molecular dynamics simulations using the Amber simulation
    package. It handles tasks such as setting up Amber-specific file paths, generating
    input files in the Amber format, constructing the command to execute Amber's
    `pmemd` engine, and managing multi-stage simulations.
    """

    @staticmethod
    def prepare_sim_config(
        simlify_config: SimlifyConfig,
    ) -> SimlifyConfig:
        r"""Preprocesses and validates the simulation configuration specifically for
        Amber simulations.

        This static method performs Amber-specific preprocessing and validation of the
        `SimlifyConfig` object. It sets up various Amber-specific file paths based on
        the configuration, such as the paths for the input file (`.in`), output file
        (`.out`), restart file (`.rst`), trajectory file (`.nc`), and MD info file
        (`.mdinfo`). It also checks for the existence of the working directory and
        creates it if necessary. For MPI-based simulations, it calculates the total
        number of CPU cores based on the SLURM configuration.

        Args:
            simlify_config: The Simlify configuration object containing options and
            parameters for the Amber simulation. This object is modified in-place.

        Returns:
            The modified Simlify configuration object with Amber-specific
                settings and paths.

        Raises:
            ValueError: If `simlify_config.label` is `None`, as it is used to construct
                file paths.
        """
        logger.info("Preparing context for Amber simulation")

        # Creates all amber-specific paths we will need.
        if simlify_config.label is None:
            raise ValueError("simlify_config.label cannot be None")
        simlify_config.engine.cli.mdin = os.path.join(
            simlify_config.run.dir_input, simlify_config.label + ".in"
        )
        simlify_config.engine.cli.mdout = os.path.join(
            simlify_config.run.dir_run, simlify_config.label + ".out"
        )
        simlify_config.engine.cli.restrt = os.path.join(
            simlify_config.run.dir_run, simlify_config.label + ".rst"
        )
        simlify_config.engine.cli.mdcrd = os.path.join(
            simlify_config.run.dir_run, simlify_config.label + ".nc"
        )
        simlify_config.engine.cli.mdinfo = os.path.join(
            simlify_config.run.dir_run, simlify_config.label + ".mdinfo"
        )

        # Checks to see if we will be writing files to a directory later.
        if not os.path.exists(simlify_config.run.dir_work):
            os.makedirs(simlify_config.dir_work, exist_ok=True)

        if "mpi" in simlify_config.engine.cli.compute_platform.lower():
            simlify_config.temp["cpu_cores"] = int(
                simlify_config.slurm.nodes * simlify_config.slurm.ntasks_per_node
            )

        return simlify_config

    @classmethod
    def get_stage_run_command(cls, simlify_config: SimlifyConfig) -> list[str]:
        r"""Prepares the bash command to run a single Amber simulation stage using the
        `pmemd` engine.

        This class method constructs the command-line invocation of Amber's `pmemd` engine
        based on the provided `SimlifyConfig`. It handles different computational platforms,
        including serial and MPI execution. If the `use_scratch` option is enabled in the
        configuration, it also adds commands to check if the simulation has already been
        completed and to move the output files from the scratch directory to the designated
        output directory.

        Args:
            simlify_config: The Simlify configuration object containing options and
                parameters for the current Amber simulation stage, such as the
                input file paths, output directory, and computational platform.

        Returns:
            A list of strings, where each string represents a line in the bash
                script that constitutes the command to run the Amber simulation stage.

        Raises:
            ValueError: If `simlify_config.label` is `None`, as it is used to construct
                file paths.

        Notes:
            This method relies on the `simlify_config` object having been preprocessed by
            the [`prepare_sim_config`][simulation.amber.prep.AmberSimPrep.prepare_sim_config]
            method.

            **Uses the following attributes from `simlify_config`:**

            - `label`: Unique label for this simulation stage.
            - `run.use_scratch`: Boolean indicating whether to use a scratch directory.
            - `run.dir_work`: Path to the main working directory.
            - `run.dir_output`: Path to the final output directory.
            - `engine.cli`: An object representing the Amber command-line interface, with
              attributes like `compute_platform`, `mdin`, `mdout`, `restrt`, and `mdcrd`.
            - `temp['cpu_cores']`: The number of CPU cores to use for MPI simulations.
        """
        if simlify_config.label is None:
            raise ValueError("label cannot be None")

        stage_commands = ["", f"echo 'Starting {simlify_config.label}'", "date"]

        if simlify_config.run.use_scratch:
            # Adds commands to check if split was already ran.
            check_path = os.path.join(
                simlify_config.run.dir_work,
                simlify_config.run.dir_output,
                simlify_config.label + ".rst",
            )
            stage_commands = ["    " + line for line in stage_commands]
            stage_commands.insert(0, f"FILE={check_path}")
            stage_commands.insert(1, 'if [ ! -f "$FILE" ]; then')

        amber_command = simlify_config.engine.cli.render()[0]
        if "mpi" in simlify_config.engine.cli.compute_platform.lower():
            amber_command = (
                f"mpirun -np {simlify_config.temp['cpu_cores']} " + amber_command
            )

        # TODO: Where do we update with the previous restart file path?
        logger.debug("Amber command: %s", amber_command[:-2])
        stage_commands.append(amber_command)

        if simlify_config.run.use_scratch:
            # Add indentation to amber command
            stage_commands[-1] = "    " + stage_commands[-1]

            # Adds commands to move scratch files to dir_output.
            stage_commands.extend(
                [
                    f"    mv {simlify_config.engine.cli.mdout} {simlify_config.run.dir_output}",
                    f"    mv {simlify_config.engine.cli.restrt} {simlify_config.run.dir_output}",
                    f"    mv {simlify_config.engine.cli.mdcrd} {simlify_config.run.dir_output}",
                    "fi",
                ]
            )

        return stage_commands

    @classmethod
    def get_stage_input_lines(cls, simlify_config: SimlifyConfig) -> list[str]:
        r"""Prepares the lines for a single Amber simulation stage's input file
        (`.in` file).

        This class method renders the content of the Amber input file based on the
        `engine.inputs` attribute of the provided `SimlifyConfig` object. The
        `engine.inputs` is expected to be an object with a `render` method that
        generates the input file content as a list of strings.

        Args:
            simlify_config (SimlifyConfig): The Simlify configuration object containing
                the input parameters for the current Amber simulation stage, typically
                defined within the `engine.inputs` attribute.

        Returns:
            list[str]: A list of strings, where each string represents a line in the
                Amber input file (`.in` file). The lines do not include the newline
                character (`\n`) at the end.
        """
        logger.info("Preparing input lines for stage {}", simlify_config.label)
        input_lines: list[str] = simlify_config.engine.inputs.render()
        return input_lines

    @classmethod
    def prepare_stage(
        cls,
        simlify_config: SimlifyConfig,
        run_commands: list[str] | None = None,
        write: bool = True,
    ) -> tuple[list[str], list[str]]:
        r"""Writes the input file for a single Amber simulation stage and builds the
        bash commands to run it.

        This class method orchestrates the preparation of a single stage in the Amber
        simulation workflow. It first retrieves the input file lines using
        `get_stage_input_lines`. If the `write` flag is True, it writes these lines to
        an Amber input file (`.in`) in the designated input directory. It then retrieves
        the bash commands to run this stage using `get_stage_run_command` and appends
        them to the cumulative list of `run_commands`. This method also handles the
        splitting of long simulations into multiple shorter segments if specified in
        the configuration.

        Args:
            simlify_config: The Simlify configuration object containing
                options and parameters for the current Amber simulation stage.
            run_commands: A list of bash commands accumulated from previous stages.
                The commands for the current stage will be appended to this list.
                Defaults to `None`, in which case a new list is initialized.
            write: A boolean flag indicating whether to write the Amber
                input file (`.in`) to disk. Defaults to `True`.

        Returns:
            A tuple containing two lists:
                -   The lines of the Amber input file generated for this stage.
                -   The updated list of bash commands, including the command(s) for this
                    stage.

        Raises:
            ValueError: If `simlify_config.label` is `None`, as it is used to construct
                file paths.

        Notes:
            This method relies on the `simlify_config` object having been preprocessed by
            the [`prepare_sim_config`][simulation.amber.prep.AmberSimPrep.prepare_sim_config]
            method.
        """
        if run_commands is None or len(run_commands) == 0:
            run_commands: list[str] = ["#!/usr/bin/env bash"]

        # We do not want to change source context in prepare_sim_config, so we do this
        # here.
        if simlify_config.run.splits > 1:
            simlify_config.engine.inputs.nstlim = int(
                simlify_config.engine.inputs.nstlim / simlify_config.run.splits
            )

        label_original = simlify_config.label
        if label_original is None:
            raise ValueError("label cannot be None")
        for i_split in range(1, simlify_config.run.splits + 1):
            if simlify_config.run.splits > 1:
                label_stage_suffix = f"_split_{i_split:03d}"  # Why n_splits < 1000
            else:
                label_stage_suffix = ""

            simlify_config.label = label_original + label_stage_suffix
            stage_input_lines = cls.get_stage_input_lines(simlify_config)
            if write:
                stage_path_input = os.path.join(
                    simlify_config.run.dir_work,
                    simlify_config.run.dir_input,
                    simlify_config.label + ".in",
                )
                logger.info("Writing input file at {}", stage_path_input)
                with open(stage_path_input, mode="w", encoding="utf-8") as f:
                    f.writelines([i + "\n" for i in stage_input_lines])

            logger.info("Adding stage's amber command to run script")
            stage_commands = cls.get_stage_run_command(simlify_config)
            run_commands.extend(stage_commands)
        simlify_config.label = label_original
        return stage_input_lines, run_commands

    @classmethod
    def prepare(cls, simlify_config: SimlifyConfig) -> None:
        """Runs all necessary steps to prepare the Amber simulations based on the
        provided configuration.

        This class method orchestrates the complete preparation process for Amber
        simulations. It handles both single-stage and multi-stage simulations as
        defined in the `simlify_config`. For multi-stage simulations, it iterates
        through the configured stages, updating the `simlify_config` with the
        parameters for each stage. For each stage, it calls `prepare_sim_config`
        to preprocess the configuration and then `prepare_stage` to generate the
        input files and run commands. Finally, if the `write` flag in the
        configuration is True, it writes the accumulated bash commands to the run
        script file specified in the configuration.

        Args:
            simlify_config: The Simlify configuration object containing
                all the necessary options and parameters for the Amber simulation workflow,
                including settings for single or multiple stages.

        Notes:
            This method assumes that the `simlify_config` object is properly initialized
            with all the required parameters for the Amber simulation.
        """
        logger.info("Prepare simulations")
        multiple_stages: bool = bool(simlify_config.stages is not None)
        if multiple_stages:
            logger.debug("There are multiple stages")
            n_stages: int = len(simlify_config.stages)
            simlify_config.update(simlify_config.stages[0])
        else:
            logger.debug("There is one stage.")
            n_stages = 1

        run_commands: list[str] = []
        for i_stage in range(1, n_stages + 1):
            logger.info("Preparing stage {}", i_stage - 1)
            simlify_config = cls.prepare_sim_config(simlify_config)
            _, run_commands = cls.prepare_stage(
                simlify_config, run_commands, simlify_config.write
            )

            simlify_config.path_restart_prev = simlify_config.path_restart
            simlify_config.path_coord_prev = simlify_config.path_coord

            if multiple_stages and i_stage < n_stages:
                simlify_config.update(simlify_config.stages[i_stage])

        if simlify_config.write:
            logger.debug("Writing run script at {}", simlify_config.path_run_write)
            with open(simlify_config.path_run_write, "w", encoding="utf-8") as f:
                f.writelines([l + "\n" for l in run_commands])
