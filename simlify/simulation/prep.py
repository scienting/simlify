from abc import ABC, abstractmethod
from collections.abc import Collection

from loguru import logger

from simlify import SimlifyConfig


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

        Notes:
            [`prepare_sim_config`][simulation.amber.prep.SimPrep.prepare_sim_config]
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
