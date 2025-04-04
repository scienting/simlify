from simlify.structure.pdb.orientation import run_minimize_box


def add_min_box_subparser(subparsers):
    r"""Command-line interface for filtering structure."""
    parser = subparsers.add_parser(
        "min-box", description="Minimize box size by rotating system."
    )
    parser.add_argument(
        "pdb_path",
        type=str,
        nargs="?",
        help="Path to PDB file",
    )
    parser.add_argument(
        "--output",
        type=str,
        nargs="?",
        help="Path to new PDB file",
    )
    parser.set_defaults(func=lambda args: cli_minimize_box(args, parser))
    return parser


def cli_minimize_box(args, parser):
    r"""Command-line interface for rotating protein to minimize box volume."""
    run_minimize_box(args.pdb_path, args.output)
