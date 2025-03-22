"""
Structure module command-line interface.
"""

from simlify.structure import extract_atoms


def add_extract_subparser(subparsers):
    r"""
    Add the `extract` subcommand to the CLI for extracting atoms or frames
    from molecular structure files.

    This subcommand wraps the `extract_atoms` function and provides arguments
    for specifying topology files, coordinate files, atom selections, frames,
    and output options.

    Args:
        subparsers: The `argparse` subparsers object to which the `extract`
            subcommand will be added.

    Returns:
        The configured `argparse` parser object for the `extract` subcommand.
    """
    parser = subparsers.add_parser(
        "extract", help="Extract atoms or frames from molecular structures"
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
        "--select",
        type=str,
        nargs="+",
        help="MDAnalysis selection string.",
    )
    parser.add_argument(
        "--frames",
        type=int,
        nargs="+",
        help="Trajectory frames.",
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
    parser.set_defaults(func=lambda args: cli_extract_atoms(args, parser))
    return parser


def cli_extract_atoms(args, parser):
    r"""
    CLI entry point for extracting atoms or frames using parsed command-line arguments.

    This function is invoked by the `extract-atoms` subcommand and acts as a
    bridge between the CLI interface and the `extract_atoms` function.

    Args:
        args: The parsed arguments from the command-line, typically from `argparse`.
        parser: The CLI parser object used to show help messages.

    Example:
        From the command-line:

        ```
        simlify structure extract mol.prmtop protein.nc \
            --select protein --frames 0 10 20 --coords traj1.nc traj2.nc --overwrite
        ```
    """
    if not any([args.topo]):
        parser.print_help()
        return

    select = args.select
    if select:
        select = " ".join(args.select)
    extract_atoms(
        path_topo=args.topo,
        select=select,
        frames=args.frames,
        path_output=args.output,
        path_coords=args.coords,
        overwrite=args.overwrite,
    )
