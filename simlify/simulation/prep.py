"""
Module defining an abstract base class for simulation preparation steps and a concrete
implementation for Amber.
"""

from typing import Any

from abc import ABC, abstractmethod
from collections.abc import Collection

from loguru import logger

from simlify.configs import SimlifyConfig


class SimPrep(ABC):
    r"""Standardized framework for preparing files for molecular simulation runs.

    This abstract base class defines the interface for preparing input files,
    run scripts, and other necessary files for executing molecular dynamics
    simulations. Concrete implementations of this class should be created for
    specific simulation packages (e.g., Amber, GROMACS) to handle their unique
    input file formats, command-line arguments, and workflow requirements.

    The `SimPrep` class provides a set of abstract methods that outline the
    essential steps involved in simulation preparation. These methods cover
    configuration preprocessing, input file generation, SLURM script creation,
    and the construction of the command to execute the simulation. By adhering
    to this framework, different simulation preparation workflows can be managed
    in a consistent and modular manner within Simlify.
    """

    def __init__(self):
        pass

    @staticmethod
    @abstractmethod
    def prepare_sim_config(
        simlify_config: SimlifyConfig,
        **kwargs: Any,
    ) -> SimlifyConfig:
        r"""Preprocesses and validates the simulation configuration for a specific
        simulation package.

        This abstract static method is responsible for performing any necessary
        preprocessing and validation of the `SimlifyConfig` object before the
        simulation preparation proceeds. This step is crucial for ensuring that all
        required parameters are set, paths are correct, and the configuration is
        consistent with the requirements of the target simulation package.
        Implementations of this method in concrete subclasses should handle the
        specific configuration needs of their respective simulation engines.

        Args:
            simlify_config: The Simlify configuration object containing
                options and parameters for the simulation. This object might be modified
                in-place during the preprocessing and validation.

        Returns:
            The potentially modified Simlify configuration object after
                preprocessing and validation.

        Raises:
            NotImplementedError: If this abstract method is called directly. Subclasses
                must provide their own implementation.
        """
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def get_stage_input_lines(
        cls,
        simlify_config: SimlifyConfig,
        **kwargs: Any,
    ) -> Collection[str]:
        r"""Prepares the lines for a single simulation stage's input file.

        This abstract class method is responsible for generating the content of the input
        file required by the simulation engine for a particular stage of the simulation.
        The generated lines should be formatted according to the input file syntax of
        the target simulation package. Implementations in concrete subclasses will define
        how the input file content is generated based on the provided `SimlifyConfig`.

        Args:
            simlify_config: The Simlify configuration object containing
                options and parameters relevant to the current simulation stage.

        Returns:
            A collection (e.g., a list) of strings, where each string represents a line
                in the input file for the simulation stage. The lines should typically
                not include the newline character (`\n`) at the end.

        Raises:
            NotImplementedError: If this abstract method is called directly. Subclasses
                must provide their own implementation.
        """
        raise NotImplementedError

    @classmethod
    def prepare_slurm(
        cls,
        file_path: str,
        simlify_config: SimlifyConfig,
        write: bool = True,
        **kwargs: Any,
    ) -> list[str]:
        r"""Prepares the lines for a SLURM submission script.

        This class method generates the content for a SLURM (Simple Linux Utility for Resource
        Management) submission script based on the SLURM configuration specified in the
        `SimlifyConfig`. The generated script can be used to submit the simulation job
        to a SLURM-managed cluster.

        Args:
            file_path: The full path to the file where the SLURM submission script
                will be written if `write: The Simlify configuration object containing
                the SLURM-specific settings (e.g., number of nodes, tasks per node,
                time limit).
            write: A boolean flag indicating whether to write the generated
                SLURM script to the specified `file_path`. Defaults to `True`.
                If `False`, the lines of the script are returned but not written to
                disk.

        Returns:
            A list of strings, where each string represents a line in the
                generated SLURM submission script, including the newline character
                (`\n`) at the end of each line.

        Notes:
            The SLURM configuration is typically defined within the
            `simlify_config.slurm` attribute, which is expected to be an object with a
            `render` method that generates the script content.
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
    def get_stage_run_command(
        cls,
        simlify_config: SimlifyConfig,
        **kwargs: Any,
    ) -> Collection[str]:
        r"""Prepares the bash command to execute a single simulation stage.

        This abstract class method is responsible for constructing the command that will
        be executed in a shell environment (e.g., bash) to run the simulation for a given
        stage. This command will typically involve invoking the simulation engine executable
        with the appropriate input file and any necessary command-line arguments. The
        specific command will vary depending on the simulation package and the configuration
        settings.

        Args:
            simlify_config: The Simlify configuration object containing
                options and parameters relevant to the current simulation stage,
                including the path to the input file and any execution settings.

        Returns:
            A collection (e.g., a list) of strings, where each string
                represents a part of the bash command to run the simulation stage. These
                parts might be joined together with spaces to form the complete command.

        Raises:
            NotImplementedError: If this abstract method is called directly. Subclasses
                must provide their own implementation.

        Notes:
            It is expected that the
            [`prepare_sim_config`][simulation.prep.SimPrep.prepare_sim_config] method
            has been called prior to this method to ensure that the `simlify_config`
            object is properly prepared.
        """
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def prepare_stage(
        cls,
        simlify_config: SimlifyConfig,
        run_commands: list[str] | None = None,
        write: bool = True,
        **kwargs: Any,
    ) -> tuple[list[str], list[str]]:
        r"""Prepares the input files for a single simulation stage and builds the
        bash commands to run it.

        This abstract class method orchestrates the preparation of a single stage in the
        simulation workflow. It typically involves generating the input file content
        using `get_stage_input_lines`, optionally writing this content to a file, and
        then constructing the command to run the simulation for this stage using
        `get_stage_run_command`.

        Args:
            simlify_config: The Simlify configuration object containing
                options and parameters relevant to the current simulation stage.
            run_commands: A list of bash commands accumulated
                from previous stages. The commands for the current stage will be
                appended to this list. Defaults to `None`, in which case a new list
                is created.
            write: A boolean flag indicating whether to write the input
                file to disk. Defaults to `True`.

        Returns:
            A tuple containing two lists:
                -   The lines of the input file generated for this stage.
                -   The updated list of bash commands, including the command(s) for this
                    stage.

        Raises:
            NotImplementedError: If this abstract method is called directly. Subclasses
                must provide their own implementation.
        """
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def prepare(cls, simlify_config: SimlifyConfig, **kwargs: Any) -> None:
        """Runs all necessary steps to prepare the simulations based on the provided configuration.

        This abstract class method serves as the main entry point for the simulation
        preparation process. It should call the other preparation methods in the correct
        order to generate all required input files, run scripts, and any other necessary
        files for the entire simulation workflow, which might consist of multiple stages.

        Args:
            simlify_config: The Simlify configuration object containing
                all the necessary options and parameters for the simulation.

        Raises:
            If this abstract method is called directly. Subclasses
                must provide their own implementation.
        """
        raise NotImplementedError
