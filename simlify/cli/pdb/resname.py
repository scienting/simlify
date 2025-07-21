# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

"""
Command-line interface for renaming residues within a Protein Data Bank (PDB) file.
"""

import argparse

from simlify.structure.pdb.names import run_replace_resnames
from simlify.structure.pdb.utils import parse_resid


def add_pdb_resname_subparser(subparsers):
    r"""Adds the 'resname' subcommand to the Simlify CLI for renaming residues.

    This function configures an `argparse` subparser named 'resname' which allows
    users to specify a PDB file, the current residue name to be replaced, the new
    residue name, and options for specifying an output file and filtering residues
    based on their indices.

    Args:
        subparsers: An `argparse._SubParsersAction` object where the 'resname'
            subparser will be added. This is typically obtained by calling
            `add_subparsers()` on an `ArgumentParser` object.

    Returns:
        The configured `argparse` parser object for the 'resname' subcommand.

    The 'resname' subcommand accepts the following arguments:

    -   **`pdb_path`:** Path to the input PDB file where residue names will be
        replaced. This argument is optional. If not provided, the
        program might expect the path through other configuration
        mechanisms.
    -   **`current_resname`:** The current three-letter code of the residue name that
        should be replaced (e.g., "ALA").
    -   **`new_resname`:** The new three-letter code to which the residue should be
        renamed (e.g., "GLY").
    -   **`--output`:** Path to the new PDB file where the modified structure
        will be saved. If not provided, the output might overwrite
        the input file or be handled in another way by the
        underlying function.
    -   **`--include`:** One or more residue indices to *only* include for renaming.
        Residues with indices not in this list will not be renamed.
        Indices should match the residue number in the PDB file.
    -   **`--exclude`:** One or more residue indices to *exclude* from renaming.
        All residues except those with indices in this list will be
        renamed. Indices should match the residue number in the
        PDB file.

    The function sets the default action for this subparser to be the
    `cli_replace_resnames` function, which will be called when the user invokes
    the 'resname' subcommand.
    """
    parser = subparsers.add_parser("resname", description="Rename residues.")
    parser.add_argument(
        "pdb_path",
        type=str,
        nargs="?",
        help="Path to PDB file",
    )
    parser.add_argument(
        "current_resname",
        type=str,
        nargs="?",
        help="Current residue name to replace",
    )
    parser.add_argument(
        "new_resname",
        type=str,
        nargs="?",
        help="New residue name",
    )
    parser.add_argument(
        "--output",
        type=str,
        nargs="?",
        help="Path to new PDB file",
    )
    parser.add_argument(
        "--include",
        type=str,
        nargs="*",
        help="Only include these residue indices.",
    )
    parser.add_argument(
        "--exclude",
        type=str,
        nargs="*",
        help="Include all residues except ones with these indices.",
    )
    parser.set_defaults(func=lambda args: cli_replace_resnames(args, parser))
    return parser


def cli_replace_resnames(
    args: argparse.Namespace, parser: argparse.ArgumentParser
) -> None:
    r"""Command-line interface function to rename residues in a PDB file.

    This function serves as the entry point when the user executes the 'resname'
    subcommand of the Simlify CLI. It receives the parsed command-line arguments
    and the argument parser object. It constructs a dictionary mapping the current
    residue name to the new residue name. It then determines whether to use a filter
    function based on the presence of the `--include` or `--exclude` arguments. If
    either of these arguments is provided, the `parse_resid` function is used as a
    filter. Finally, it calls the `run_replace_resnames` function from the
    `simlify.structure.pdb.names` module to perform the residue renaming operation.

    Args:
        args: An `argparse.Namespace` object containing the parsed command-line
            arguments for the 'resname' subcommand.
        parser: The argument parser object for the
            'resname' subcommand, used to display help messages if necessary (though
            not explicitly used in this function).

    The function retrieves the following arguments from the `args` object:

    - `pdb_path`: Path to the input PDB file.
    - `current_resname`: The residue name to be replaced.
    - `new_resname`: The new residue name.
    - `output`: Path to the output PDB file.
    - `include`: A list of residue indices to include for renaming.
    - `exclude`: A list of residue indices to exclude from renaming.

    It creates a `resname_map` dictionary with the provided current and new
    residue names. It also sets the `fn_filter` to `parse_resid` if either
    `args.include` or `args.exclude` is provided, otherwise it remains `None`.
    Finally, it calls [`run_replace_resnames`][structure.pdb.names.run_replace_resnames]
    with these parameters to perform the residue renaming.
    """
    resname_map = {args.current_resname: args.new_resname}
    if (args.include is not None) or (args.exclude is not None):
        fn_filter = parse_resid
    else:
        fn_filter = None
    run_replace_resnames(
        args.pdb_path, resname_map, args.output, fn_filter, args.include, args.exclude
    )
