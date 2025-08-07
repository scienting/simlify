# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

"""
Command-line interface for unifying residue and atom names of water molecules in a PDB file.
"""

import argparse

from simlify.structure.pdb.names import run_unify_water_labels


def add_pdb_water_subparser(subparsers):
    r"""Adds the 'unify-wat' subcommand to the Simlify CLI for unifying water labels.

    This function configures an `argparse` subparser named 'unify-wat' which allows
    users to specify a PDB file and an optional output path. The 'unify-wat'
    command will standardize the residue and atom names of all water molecules
    within the PDB file to a consistent format (e.g., residue name to 'HOH' and
    oxygen atom name to 'O').

    Args:
        subparsers: An `argparse._SubParsersAction` object where the 'unify-wat'
            subparser will be added. This is typically obtained by calling
            `add_subparsers()` on an `ArgumentParser` object.

    Returns:
        The configured `argparse` parser object for the 'unify-wat' subcommand.

    The 'unify-wat' subcommand accepts the following arguments:

    -   **`pdb_path`:** Path to the input PDB file containing water molecules whose
        labels will be unified. This argument is optional. If not
        provided, the program might expect the path through other
        configuration mechanisms.
    -   **`--output:** Path to the new PDB file where the water labels have been
        unified. If not provided, the output might overwrite the
        input file or be handled in another way by the underlying
        function.

    The function sets the default action for this subparser to be the
    `cli_unify_water_labels` function, which will be called when the user invokes
    the 'unify-wat' subcommand.
    """
    parser = subparsers.add_parser(
        "unify-wat", description="Unify water residue and atom names"
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
    parser.set_defaults(func=lambda args: cli_unify_water_labels(args, parser))
    return parser


def cli_unify_water_labels(
    args: argparse.Namespace, parser: argparse.ArgumentParser
) -> None:
    r"""Command-line interface function to unify water residue and atom names in a PDB file.

    This function serves as the entry point when the user executes the 'unify-wat'
    subcommand of the Simlify CLI. It receives the parsed command-line arguments
    and the argument parser object. It extracts the paths to the input PDB file
    and the desired output PDB file and then calls the `run_unify_water_labels`
    function from the `simlify.structure.pdb.names` module to perform the
    unification of water labels.

    Args:
        args: An `argparse.Namespace` object containing the parsed command-line
            arguments for the 'unify-wat' subcommand.
        parser (argparse.ArgumentParser): The argument parser object for the
            'unify-wat' subcommand. Although passed, it's not directly used in
            this specific implementation.

    Returns:
        None

    The function retrieves the following arguments from the `args` object:
        - `pdb_path`: Path to the input PDB file.
        - `output`: Path to the output PDB file.

    It then calls the
    [`run_unify_water_labels`][structure.pdb.names.run_unify_water_labels]
    function with these paths to perform the standardization of water residue and
    atom names.
    """
    run_unify_water_labels(args.pdb_path, output_path=args.output)
