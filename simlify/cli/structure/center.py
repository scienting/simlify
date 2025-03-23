from simlify.structure.center import run_center_structure


def add_center_subparser(subparsers):
    r"""Command-line interface for centering structure."""
    parser = subparsers.add_parser(
        "center", description="Move center of mass of all frames to the origin."
    )
    parser.add_argument(
        "topo",
        type=str,
        nargs="?",
        help="Path to topology file.",
    )
    parser.add_argument(
        "output",
        type=str,
        nargs="?",
        help="Path to coordinate file to save.",
    )
    parser.add_argument(
        "--coords",
        type=str,
        nargs="*",
        help="Paths to coordinate files. Will concatenate these into one file.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite coordinate file.",
    )
    parser.set_defaults(func=lambda args: cli_center_atoms(args, parser))
    return parser


def cli_center_atoms(args, parser):
    if not any([args.topo]):
        parser.print_help()
        return

    select = args.select
    if select:
        select = " ".join(args.select)
    run_center_structure(
        path_topo=args.topo,
        path_output=args.output,
        path_coords=args.coords,
        overwrite=args.overwrite,
    )
