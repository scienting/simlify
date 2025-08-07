# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

"""
Command-line interface for unifying numbering within a Protein Data Bank (PDB) file.
"""

import argparse

from simlify.structure.pdb.numbering.main import run_unify_numbering


def add_pdb_unify_numbering_subparser(subparsers):
    r"""Adds the `unify_numbering` subcommand to the Simlify CLI for unifying atom
    and residue numbering in PDB files.

    This function configures an `argparse` subsparser named `unify_numbering` which
    allows users to specify a PDB file and option to control whether the initial
    residue ID should be reset to 1.

    Args:
        subparsers: An `argparse._SubParsersAction` object where the `unify_numbering`
            subparser will be added. This is typically obtained by calling
            `add_subparsers()` on an `ArgumentParser` object.
    Returns:
        argparse.ArgumentParser: The configured `argparse` parser object for the
        `unify_numbering` subcommand.

    The `unify_numbering` subcommand accepts the following arguments:

    positional arguments:
        pdb_path    Path to the input PDB file whose numbering will be processed.
                    This argument is optional. If not provided, the program might
                    expect the path through other configuration mechanisms.

    optional arguments:
        --output            Path to the new PDB file where the output will be saved.
                            If not provided, the output might overwrite the input file
                            or be handled in another way by the underlying function.
        --reset_init_resid   reset original first residue ID; do not reset to 1. This
                            Not modify the numbering in the initial amino acid in the
                            chain.

    The function sets the default action for this subparser to be the
    `cli_unify_numbering` function, which will be called when the user invokes
    the `unify_numbering` subcommand.
    """

    parser = subparsers.add_parser(
        "unify_numbering", description="Unify atom and residue IDs in PDB"
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
    parser.add_argument(
        "--reset_init_resid",
        action="store_true",
        help="Do not reset first residue ID to 1.",
    )
    parser.set_defaults(func=lambda args: cli_unify_numbering(args, parser))
    return parser


def cli_unify_numbering(
    args: argparse.Namespace, parser: argparse.ArgumentParser
) -> None:
    r"""Command-line interface function to for unifying atom and residue numbering in
    PDB files.

    This function serves as the entry point when the user executes the `unify_numbering`
    subcommand of the Simlify CLI. It receives the parsed command-line arguments
    and the argument parser object. It extracts the paths to the input PDB file
    and the desired output PDB file and then calls the `run_unify_numbering` function
    from the `simlify.structure.pdb.numbering.main` to renumber resIDs and atom numbers.
    Args:
        args: An `argparse.Namespace` object containing the parsed command-line arguments for
            the 'unify_numbering' subcommand.
        parser (argparse.ArgumentParser): The argument parser object for the
            `unify_numbering` subcommand, used to display help messages if necessary (though
            not explicitly used in this function).

    Returns:
        None

    The function retrieves the following arguments from the `args` object:
        - `pdb_path` is the path to the PDB file to be processed.
        - `--reset_init_resid` is an optional flag. If present, the initial residue
        ID will not be reset to 1.

    It then directly calls the `run_unify_numbering` function with these paths to
    perform the renumbering operation.
    See Also:
        `run_unify_numbering`: The underlying function that performs the PDB numbering
        unification.
    """
    if args.pdb_path is None:
        raise RuntimeError("Input must be specified")
    run_unify_numbering(args.pdb_path, args.output, args.reset_init_resid)
