from pydantic import BaseModel

from simlify.schemas.utils import YamlIO


class RunConfig(BaseModel, YamlIO):
    """Configuration options during the simulation runtime."""

    dir_work: str = "."
    """
    Working directory during runtime. This can be the current work directory or
    a workload manager directory such as `$SLURM_SUBMIT_DIR`.
    """

    dir_run: str = "."
    """
    Directory that calculations are performed in. This can be the same as
    [`dir_work`][schemas.RunConfig.dir_work]
    or some scratch space like `$SLURM_SCRATCH`.
    """

    dir_input: str = "."
    """
    Path to a directory relative to
    [`dir_work`][schemas.RunConfig.dir_work] that will contain
    input files.
    """

    dir_output: str = "."
    """
    Path to a directory relative to
    [`dir_work`][schemas.RunConfig.dir_work] that the simulation will
    store output files.
    """

    splits: int = 1
    """Split simulation stage into several chunks."""

    use_scratch: bool = False
    """
    Treat [`dir_run`][schemas.RunConfig.dir_run] as a scratch directory
    by copying input files to that location, running the simulation there, and
    then copying output files back.
    """
