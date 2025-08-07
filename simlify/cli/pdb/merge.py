# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

"""
Command-line interface for merging multiple Protein Data Bank (PDB) files into a single file.
"""

import argparse

from simlify.structure.pdb.utils import run_merge_pdbs


def add_pdb_merge_subparser(subparsers):
    r"""Adds the 'merge' subcommand to the Simlify CLI for combining multiple PDB files.

    This function configures an `argparse` subparser named 'merge' which allows
    users to specify one or more PDB files to be merged together. It also provides
    an option to specify the output path for the merged PDB file.

    Args:
        subparsers: An `argparse._SubParsersAction` object where the 'merge'
            subparser will be added. This is typically obtained by calling
            `add_subparsers()` on an `ArgumentParser` object.

    Returns:
        The configured `argparse` parser object for the 'merge' subcommand.

    The 'merge' subcommand accepts the following arguments:

    -   **`pdb_paths`:** One or more paths to the PDB files that should be merged.
        These files will be concatenated in the order they are provided.
    -   **`--output`:**  Path to the new PDB file where the content of all input
        PDB files will be written.

    The function sets the default action for this subparser to be the
    `cli_merge_pdbs` function, which will be called when the user invokes
    the 'merge' subcommand.
    """
    parser = subparsers.add_parser("merge", description="Merge PDB files")
    parser.add_argument(
        "pdb_paths",
        type=str,
        nargs="+",
        help="PDB files to merge",
    )
    parser.add_argument(
        "--output",
        type=str,
        nargs="?",
        help="Path to new PDB file",
    )
    parser.set_defaults(func=lambda args: cli_merge_pdbs(args, parser))
    return parser


def cli_merge_pdbs(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    r"""Command-line interface function to merge multiple PDB files into one.

    This function serves as the entry point when the user executes the 'merge'
    subcommand of the Simlify CLI. It receives the parsed command-line arguments
    and the argument parser object. It first checks if the `--output` argument
    was provided, as it is a requirement for this operation. If not, it raises
    a `RuntimeError`. Otherwise, it extracts the list of input PDB file paths
    and the output path and calls the `run_merge_pdbs` function from the
    `simlify.structure.pdb.utils` module to perform the merging operation.

    Args:
        args: An `argparse.Namespace` object containing the parsed command-line
            arguments for the 'merge' subcommand.
        parser: The argument parser object for the 'merge' subcommand, used to display
            help messages if necessary (though not explicitly used in this function).

    Raises:
        RuntimeError: If the `--output` argument is not specified by the user.

    The function retrieves the following arguments from the `args` object:

    - `pdb_paths`: A list of paths to the PDB files to be merged.
    - `output`: The path to the output PDB file.

    It then unpacks the `pdb_paths` list as positional arguments to the
    [`run_merge_pdbs`][structure.pdb.utils.run_merge_pdbs] function and provides
    the `output` path as a keyword argument.
    """
    if args.output is None:
        raise RuntimeError("--output must be specified")
    run_merge_pdbs(*args.pdb_paths, output_path=args.output)
