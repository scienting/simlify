from simlify import SimlifyConfig
from simlify.simulation.slurm import run_slurm_prep


def add_slurm_subparser(subparsers):
    parser = subparsers.add_parser("slurm", description="Render and write slurm files.")
    parser.add_argument(
        "dir_write",
        type=str,
        nargs="?",
        help="Directory to write input files.",
    )
    parser.add_argument(
        "path_slurm_write",
        type=str,
        nargs="?",
        help="Local path to write a slurm submission script.",
    )
    parser.add_argument(
        "--yaml",
        type=str,
        nargs="*",
        help="Paths to YAML files to use in decreasing precedence.",
    )
    parser.set_defaults(func=lambda args: cli_extract_atoms(args, parser))
    return parser


def cli_extract_atoms(args, parser):
    simlify_config = SimlifyConfig()
    if args.yaml is None:
        args.yaml = []
    for yaml_path in reversed(args.yaml):
        simlify_config.from_yaml(yaml_path)
    simlify_config.slurm.job_name = args.name_job
    run_slurm_prep(
        args.dir_write,
        args.path_slurm_write,
        simlify_config,
    )
