from simlify import SimlifyConfig
from simlify.simulation.topo import run_gen_topo


def add_topo_subparser(subparsers):
    parser = subparsers.add_parser("topo", description="Prepare simulation topology.")
    parser.add_argument(
        "structure",
        type=str,
        nargs="?",
        help="Path to a structure file.",
    )
    parser.add_argument(
        "import_string",
        type=str,
        nargs="?",
        help="Import string to a topology generation class.",
    )
    parser.add_argument(
        "--yaml",
        type=str,
        nargs="+",
        help="Paths to YAML files to use in decreasing precedence.",
    )
    parser.add_argument(
        "--work",
        type=str,
        nargs="?",
        help="Work directory for preparing calculations. This overrides the YAML files.",
        default=None,
    )
    parser.set_defaults(func=lambda args: cli_topo(args, parser))
    return parser


def cli_topo(args, parser):
    simlify_config = SimlifyConfig()
    if args.yaml is None:
        args.yaml = []
    for yaml_path in reversed(args.yaml):
        simlify_config.from_yaml(yaml_path)
    if args.work is not None:
        simlify_config.run.dir_work = args.work
    run_gen_topo(
        args.structure,
        args.import_string,
        simlify_config,
    )
