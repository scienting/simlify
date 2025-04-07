"""Defines the command-line interface for centering molecular structures."""

import argparse

from simlify.structure.center import run_center_structure


def add_center_subparser(subparsers):
    r"""Adds the 'center' subcommand to the provided subparsers.

    This function creates a subparser named 'center' under the main Simlify CLI,
    allowing users to move the center of mass of all frames within a molecular
    structure to the origin (0.0, 0.0, 0.0). It defines the necessary command-line
    arguments for specifying the topology, coordinate files, output path, and
    overwrite behavior.

    Args:
        subparsers: An `argparse._SubParsersAction` object where the 'center'
            subparser will be added. This is typically obtained by calling
            `add_subparsers()` on an `ArgumentParser` object.

    Returns:
        argparse.ArgumentParser: The newly created 'center' subparser.

    The 'center' subparser accepts the following arguments:

    positional arguments:
        topo        Path to topology file.
        output      Path to coordinate file to save the centered structure.

    optional arguments:
        --coords    Paths to one or more coordinate files. If multiple files are
                    provided, they will be concatenated and treated as frames of
                    the same trajectory.
        --overwrite Overwrite the output coordinate file if it already exists.

    The function sets the default action for this subparser to be the
    `cli_center_atoms` function, which will be called when the user invokes
    the 'center' subcommand.
    """
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


def cli_center_atoms(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    """Handles the command-line invocation of the structure centering functionality.

    This function serves as the entry point when the user executes the 'center'
    subcommand of the Simlify CLI. It takes the parsed command-line arguments
    and the argument parser object as input. It performs a basic check to ensure
    that at least a topology file is provided. If not, it prints the help message
    for the 'center' subcommand and returns. Otherwise, it extracts the relevant
    arguments and calls the underlying `run_center_structure` function to perform
    the structure centering operation.

    Args:
        args (argparse.Namespace): An object containing the parsed command-line
            arguments for the 'center' subcommand.
        parser (argparse.ArgumentParser): The argument parser object for the
            'center' subcommand, used to display help messages if necessary.

    Returns:
        None

    The function retrieves the following arguments from the `args` object:
        - `topo`: Path to the topology file.
        - `output`: Path to the output coordinate file.
        - `coords`: A list of paths to coordinate files.
        - `overwrite`: A boolean indicating whether to overwrite the output file.

    It then calls the `run_center_structure` function with these arguments to
    perform the actual centering of the molecular structure(s).
    """
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
