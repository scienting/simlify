# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

"""
Command-line interface for extracting atoms or frames from molecular structure files.
"""

import argparse

from simlify.structure import extract_atoms


def add_extract_subparser(subparsers):
    r"""
    Adds the `extract` subcommand to the Simlify CLI for extracting atoms or frames
    from molecular structure files.

    This function configures an `argparse` subparser named 'extract' which allows
    users to specify the topology and coordinate files, define atom selections using
    MDAnalysis selection syntax, specify trajectory frames to extract, and set the
    output file path and overwrite options.

    Args:
        subparsers: An `argparse._SubParsersAction` object where the `extract`
            subcommand will be added. This is typically obtained by calling
            `add_subparsers()` on an `ArgumentParser` object.

    Returns:
        argparse.ArgumentParser: The configured `argparse` parser object for the
        `extract` subcommand.

    The `extract` subcommand accepts the following arguments:

    positional arguments:
        topo        Path to the topology file (e.g., PDB, PRMTOP). This argument
                    is optional, but at least one of `topo` or `--coords` must be
                    provided.
        output      Path to the coordinate file where the extracted atoms or
                    frames will be saved (e.g., output.pdb, output.nc).

    optional arguments:
        --select    One or more strings specifying the atom selection using the
                    MDAnalysis selection language. Multiple strings will be joined
                    by spaces. This option allows you to select specific atoms
                    based on various criteria (e.g., "protein", "resid 1-10",
                    "name CA").
        --frames    One or more integers specifying the indices of the trajectory
                    frames to extract. If this option is not provided, all frames
                    will be processed.
        --coords    One or more paths to coordinate files (e.g., trajectory files
                    like TRR, DCD, XTC, NetCDF). If multiple coordinate files are
                    provided, they will be concatenated and treated as a single
                    trajectory. At least one coordinate file must be provided if
                    no topology file is given.
        --overwrite If this flag is present, the output coordinate file will be
                    overwritten if it already exists. By default, the program will
                    prevent overwriting existing files.

    The function sets the default action for this subparser to be the
    `cli_extract_atoms` function, which will be called when the user invokes
    the 'extract' subcommand.
    """
    parser = subparsers.add_parser(
        "extract", description="Extract atoms or frames from molecular structures"
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


def cli_extract_atoms(
    args: argparse.Namespace, parser: argparse.ArgumentParser
) -> None:
    r"""
    Command-line interface function to extract atoms or frames from molecular structures.

    This function serves as the entry point when the user invokes the `extract`
    subcommand of the Simlify CLI. It receives the parsed command-line arguments
    and the argument parser object. It performs a basic validation to ensure that
    either a topology file (`topo`) or at least one coordinate file (`coords`) is
    provided. If neither is present, it prints the help message for the `extract`
    subcommand and exits. Otherwise, it prepares the arguments and calls the
    `extract_atoms` function from the `simlify.structure` module to perform the
    actual extraction.

    Args:
        args: An `argparse.Namespace` object containing the parsed command-line
            arguments for the `extract` subcommand.
        parser: The `argparse.ArgumentParser` object for the `extract` subcommand,
            used to display help messages if necessary.

    Returns:
        None

    The function retrieves the following arguments from the `args` object:
        - `topo`: Path to the topology file.
        - `output`: Path to the output coordinate file.
        - `select`: A list of strings representing the MDAnalysis atom selection.
        - `frames`: A list of integers representing the trajectory frames to extract.
        - `coords`: A list of paths to coordinate files.
        - `overwrite`: A boolean indicating whether to overwrite the output file.

    It then processes the `select` argument by joining the list of strings into a
    single selection string. Finally, it calls the `extract_atoms` function with
    the extracted and processed arguments to perform the atom or frame extraction.

    Example Usage:
        To extract protein atoms from frames 0, 10, and 20 of a trajectory:

        ```bash
        simlify structure extract topology.prmtop output.pdb \
            --select protein --frames 0 10 20 --coords trajectory.nc
        ```

        To extract all atoms from a specific residue range and overwrite the output file:

        ```bash
        simlify structure extract input.pdb extracted.pdb \
            --select 'resid 1 to 50' --overwrite
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
