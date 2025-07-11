from pydantic import BaseModel

from simlify.configs.utils import Render, YamlIO


class SlurmConfig(BaseModel, YamlIO, Render):
    """Context manager for Slurm job submission scripts.

    This class provides a structured way to define and manage the configuration for
    submitting jobs to a Slurm workload manager. Each attribute corresponds to a
    specific Slurm configuration parameter or job setup step.
    """

    def render(self, with_newlines: bool = False) -> list[str]:
        lines: list[str] = ["#!/bin/bash", ""]

        for key, value in self.model_dump(exclude_none=True).items():
            # Skip keys that we have to manually handle
            if key in (
                "modules",
                "env_vars",
                "commands_pre",
                "commands_run",
                "commands_post",
            ):
                continue

            key = key.replace("_", "-")
            line = f"#SBATCH --{key}={value}"
            lines.append(line)

        lines.extend(self.commands_pre)
        for module in self.modules:
            lines.append(f"module load {module}")
        for key, var in self.env_vars.items():
            lines.append(f"export {key}={var}")
        lines.extend(self.commands_run)
        lines.extend(self.commands_post)
        lines.append("")

        if with_newlines:
            lines = [line + "\n" for line in lines]

        return lines

    job_name: str = "job"
    """Specify a name for the job allocation. The specified name will appear along
    with the job id number when querying running jobs on the system. The default is
    the name of the batch script, or just "sbatch" if the script is read on sbatch's
    standard input.

    [More information](https://slurm.schedmd.com/sbatch.html#OPT_job-name)

    Example:
        `"data_analysis_job"`
    """

    nodes: int = 1
    """The minimum number of nodes to use for the Slurm job.

    Adjust this based on the job's resource requirements. For instance, a large parallel
    job might need several nodes, while a smaller job might only need one.

    [More information](https://slurm.schedmd.com/sbatch.html#OPT_nodes)

    Example:
        `4`
    """

    ntasks_per_node: int | None = None
    """Request that `ntasks` be invoked on each node.

    This typically corresponds to the number of CPU cores to use on each node. Adjust
    this based on the node's capabilities and the parallelism of your job.

    [More information](https://slurm.schedmd.com/sbatch.html#OPT_ntasks-per-node)

    Example:
        `16`
    """

    ntasks: int | None = None
    """sbatch does not launch tasks, it requests an allocation of resources and
    submits a batch script. This option advises the Slurm controller that job steps
    run within the allocation will launch a maximum of number tasks and to provide
    for sufficient resources.

    [More information](https://slurm.schedmd.com/sbatch.html#OPT_ntasks)

    Example:
        `16`
    """

    output: str = "slurm-%j.out"
    """Path for the job's standard output file.

    Use `%j` to include the job ID in the filename, ensuring that each job's output is
    saved to a unique file.

    Example:
        `"logs/job_output_%j.out"`
    """

    error: str = "slurm-%j.err"
    """Path for the job's error output file.

    Similar to `output_path`, using `%j` in the filename ensures that errors for each
    job are logged separately.

    Example:
        `"logs/job_errors_%j.err"`
    """

    time: str = "1-00:00:00"
    """Maximum time for the job.

    Specified in the format `D-HH:MM:SS`. Adjust this based on the expected runtime of
    your job.

    Example:
        `"0-12:00:00"` for a 12-hour job.
    """

    cluster: str = "smp"
    """Cluster name where the job will run.

    Ensure this matches the available cluster names in your Slurm environment. This
    helps direct the job to the appropriate set of resources.

    Example:
        `"hpc_cluster"`
    """

    partition: str = "smp"
    """Partition name to submit the job to.

    Choose an appropriate partition based on resource needs and availability. Partitions
    can have different resource limits and policies.

    [More information](https://slurm.schedmd.com/sbatch.html#OPT_partition)

    Example:
        `"short"`
    """

    account: str | None = None
    """Charge resources used by this job to specified account. The account is an
    arbitrary string.

    Set this to your project's account name to properly attribute resource usage, or
    leave it as `None` if not needed.

    [More information](https://slurm.schedmd.com/sbatch.html#OPT_account)

    Example:
        `"research_project_123"`
    """

    constraint: str | None = None
    """Nodes can have features assigned to them by the Slurm administrator. Users can
    specify which of these features are required by their job using the constraint
    option.

    [More information](https://slurm.schedmd.com/sbatch.html#OPT_constraint)
    """

    cores_per_socket: int | None = None
    """Restrict node selection to nodes with at least the specified number of cores per
    socket.

    [More information](https://slurm.schedmd.com/sbatch.html#OPT_cores-per-socket)
    """

    cpus_per_gpu: int | None = None
    """Request that ncpus processors be allocated per allocated GPU. Steps inheriting
    this value will imply `--exact`. Not compatible with the
    [`--cpus-per-task`][configs.workflow.slurm.SlurmConfig.cpus_per_task] option.

    [More information](https://slurm.schedmd.com/sbatch.html#OPT_cpus-per-gpu)
    """

    cpus_per_task: int | None = None
    """Advise the Slurm controller that ensuing job steps will require `ncpus` number
    of processors per task. Without this option, the controller will just try to
    allocate one processor per task.

    [More information](https://slurm.schedmd.com/sbatch.html#OPT_cpus-per-task)
    """

    gpus: int | str | None = None
    """Specify the total number of GPUs required for the job.

    [More information](https://slurm.schedmd.com/sbatch.html#OPT_gpus)
    """

    gpus_per_node: int | str | None = None
    """Specify the number of GPUs required for the job on each node included in the
    job's resource allocation.

    [More information](https://slurm.schedmd.com/sbatch.html#OPT_gpus-per-node)
    """

    gpus_per_socket: int | str | None = None
    """Specify the number of GPUs required for the job on each socket included in the
    job's resource allocation.

    [More information](https://slurm.schedmd.com/sbatch.html#OPT_gpus-per-socket)
    """

    gpus_per_task: int | str | None = None
    """Specify the number of GPUs required for the job on each task to be spawned in
    the job's resource allocation.

    [More information](https://slurm.schedmd.com/sbatch.html#OPT_gpus-per-task)
    """

    gres: str | None = None
    """Specifies a comma-delimited list of generic consumable resources. The format
    for each entry in the list is `"name[[:type]:count]"`. The name is the type
    of consumable resource (e.g. gpu). The type is an optional classification for
    the resource (e.g. a100). The count is the number of those resources with a
    default value of 1.

    [More information](https://slurm.schedmd.com/sbatch.html#OPT_gres)
    """

    mem: str | None = None
    """Specify the real memory required per node. Default units are megabytes.
    Different units can be specified using the suffix [K|M|G|T].

    [More information](https://slurm.schedmd.com/sbatch.html#OPT_mem)
    """

    mem_per_cpu: str | None = None
    """Minimum memory required per usable allocated CPU. Default units are megabytes.
    The default value is `DefMemPerCPU` and the maximum value is `MaxMemPerCPU`.

    [More information](https://slurm.schedmd.com/sbatch.html#OPT_mem-per-cpu)
    """

    modules: list[str] = []
    """List of modules to load before running the job.

    Include all necessary software modules that your job requires. This ensures the
    environment is correctly set up before execution.

    Example:
        `["python/3.8", "gcc/9.2"]`
    """

    env_vars: dict[str, str] = {}
    """Dictionary of environment variables to set before running the job.

    Use this to configure the job's environment, setting any necessary environment
    variables.

    Example:
        `{"OMP_NUM_THREADS": "16", "MY_VARIABLE": "value"}`
    """

    commands_pre: list[str] = []
    """List of commands to run before the main job command.

    Useful for setup tasks like copying files, creating directories, or loading
    additional software. These commands will be executed before the main job starts.

    Example:
        `["mkdir -p /scratch/my_job", "cp input.dat /scratch/my_job/"]`
    """

    commands_run: list[str] = []
    """List of main commands to run for the job.

    This should include the primary executable or script for the job. These are the
    main tasks that the job will perform.

    Example:
        `["python my_script.py", "./run_simulation.sh"]`
    """

    commands_post: list[str] = []
    """List of commands to run after the main job command.

    Use this for cleanup tasks or additional processing. These commands will be executed
    after the main job tasks are completed.

    Example:
        `["cp /scratch/my_job/output.dat .", "rm -rf /scratch/my_job"]`
    """
