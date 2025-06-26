"""
Command-line interface for filtering specific records within a PDB file.
"""

import argparse

from loguru import logger

from simlify.structure.pdb.utils import run_filter_pdb


def add_pdb_filter_subparser(subparsers):
    r"""Adds the 'filter' subcommand to the Simlify CLI for filtering PDB file lines.

    This function configures an `argparse` subparser named 'filter' which allows
    users to specify a PDB file, an optional output path, and a list of PDB record
    types to keep in the output file. Lines with record types not specified will
    be removed.

    Args:
        subparsers: An `argparse._SubParsersAction` object where the 'filter'
            subparser will be added. This is typically obtained by calling
            `add_subparsers()` on an `ArgumentParser` object.

    Returns:
        The configured `argparse` parser object for the 'filter' subcommand.

    The 'filter' subcommand accepts the following arguments:

    -   **`pdb_path`:** Path to the input PDB file that will be filtered. This
        argument is optional. If not provided, the program might
        expect the path through other configuration mechanisms.
    -   **`--output`:** Path to the new PDB file where the filtered content will
        be saved. If not provided, the output might overwrite the
        input file or be handled in another way by the underlying
        function.
    -   **`--records`:** One or more PDB record types (e.g., ATOM, HETATM, TER) to
        keep in the output file. Only lines starting with these
        record types will be included in the output. If no records
        are specified, the behavior will depend on the underlying
        `run_filter_pdb` function.

    The function sets the default action for this subparser to be the
    `cli_filter_pdb` function, which will be called when the user invokes
    the 'filter' subcommand.
    """
    parser = subparsers.add_parser("filter", description="Filter PDB lines")
    parser.add_argument(
        "pdb_path",
        type=str,
        nargs="?",
        help="Path to PDB file",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        nargs="?",
        help="Path to new PDB file",
    )
    parser.add_argument(
        "--records",
        type=str,
        help="PDB records to keep. For example: ATOM,HETATM,MODEL",
    )
    parser.set_defaults(func=lambda args: cli_filter_pdb(args, parser))
    return parser


def cli_filter_pdb(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    r"""Command-line interface function to filter lines in a PDB file.

    This function serves as the entry point when the user executes the 'filter'
    subcommand of the Simlify CLI. It receives the parsed command-line arguments
    and the argument parser object. It extracts the paths to the input PDB file
    and the desired output PDB file, as well as the list of record types to keep.
    It then calls the `run_filter_pdb` function from the
    `simlify.structure.pdb.utils` module to perform the filtering operation.

    Args:
        args: An `argparse.Namespace` object containing the parsed command-line
            arguments for the 'filter' subcommand.
        parser: The argument parser object for the
            'filter' subcommand, used to display help messages if necessary (though
            not explicitly used in this function).

    The function retrieves the following arguments from the `args` object:

    - `pdb_path`: Path to the input PDB file.
    - `output`: Path to the output PDB file.
    - `records`: A list of PDB record types to keep.

    It then directly calls the [`run_filter_pdb`][structure.pdb.utils.run_filter_pdb]
    function with these arguments to perform the filtering of the PDB file.
    """
    records = args.records
    if records is not None:
        records = tuple(r.strip() for r in args.records.split(","))
    pdb_filtered = run_filter_pdb(args.pdb_path, args.output, records)
    if args.output is None:
        print("".join(pdb_filtered))
