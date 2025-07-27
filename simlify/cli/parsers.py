# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

"""Defines the argument parsing structure for the Simlify command-line interface (CLI).

This module utilizes the `argparse` library to create a hierarchical argument parsing
system for the Simlify CLI. It defines the top-level arguments and sets up subparsers
for different categories of functionalities, such as structure manipulation, PDB file
operations, and simulation preparation. Each subparser further organizes commands
related to its specific domain. This modular approach allows for a well-structured
and user-friendly command-line interface for Simlify.
"""

import argparse
import sys

from simlify.cli.pdb.filter import add_pdb_filter_subparser
from simlify.cli.pdb.merge import add_pdb_merge_subparser
from simlify.cli.pdb.numbering import add_pdb_unify_numbering_subparser
from simlify.cli.pdb.resname import add_pdb_resname_subparser
from simlify.cli.pdb.water import add_pdb_water_subparser
from simlify.cli.prep.slurm import add_slurm_subparser
from simlify.cli.prep.topo import add_topo_subparser
from simlify.cli.structure.center import add_center_subparser
from simlify.cli.structure.extract import add_extract_subparser
from simlify.cli.structure.min_box import add_min_box_subparser


def create_parser() -> argparse.ArgumentParser:
    """Creates and configures the main argument parser for the Simlify CLI.

    This function initializes an `ArgumentParser` object with a description of the
    Simlify CLI. It then defines the top-level arguments that are common across
    various subcommands, such as verbosity flags for logging, a log file option,
    and a configuration file path. Additionally, it sets up subparsers to handle
    different categories of commands, including 'structure', 'pdb', and 'prep',
    each with their own nested subcommands.

    Returns:
        The configured argument parser for the Simlify CLI.

    Notes:
        The function defines the following top-level arguments:

        -   `-v`: Enables more verbose logging.
        -   `-vv`: Enables even more verbose logging.
        -   `--logfile`: Specifies a file to write log messages to.
        -   `--config`: Specifies the path to a YAML configuration file to load settings
            from.

        It also creates subparsers for the following top-level commands:

        -   `structure`: Contains commands for manipulating molecular structure files,
            including:
            - `extract`: Subcommand for extracting parts of a structure.
            - `center`: Subcommand for centering a structure.
            - `min_box`: Subcommand for adjusting the size of the simulation box.
        -   `pdb`: Contains commands for operating on Protein Data Bank (PDB) files,
            including:
            - `filter`: Subcommand for filtering lines in a PDB file.
            - `merge`: Subcommand for merging multiple PDB files.
            - `resname`: Subcommand for modifying residue names in a PDB file.
            - `water`: Subcommand for manipulating water molecules in a PDB file.
            - `unify_numbering` : Subcommand for changing residue and atomic numbers.
        -   `prep`: Contains commands for preparing simulation files, including:
            - `slurm`: Subcommand for preparing SLURM submission scripts.
            - `topo`: Subcommand for generating topology files.

        The function utilizes helper functions (e.g., `add_extract_subparser`) from other
        modules within the `simlify.cli` package to define the arguments and actions
        for each specific subcommand.
    """
    parser = argparse.ArgumentParser(description="Simlify CLI")
    parser.add_argument("-v", action="store_true", help="More log verbosity.")
    parser.add_argument("-vv", action="store_true", help="Even more log verbosity.")
    parser.add_argument("--logfile", help="Specify a file to write logs to.")
    parser.add_argument("--config", help="Path to YAML configuration file.")

    # Top-level subparsers
    subparsers = parser.add_subparsers(dest="command", help="Top-level commands")

    # Structure subcommand and its nested commands
    structure_parser = subparsers.add_parser(
        "structure", help="Commands for molecular structure operations"
    )
    parser._structure_parser = structure_parser  # type: ignore
    structure_subparsers = structure_parser.add_subparsers(
        dest="structure_command", help="Structure commands"
    )
    add_extract_subparser(structure_subparsers)
    add_center_subparser(structure_subparsers)
    add_min_box_subparser(structure_subparsers)

    # PDB subcommand and its nested commands
    pdb_parser = subparsers.add_parser("pdb", help="Commands for PDB file operations")
    parser._pdb_parser = pdb_parser  # type: ignore
    pdb_subparsers = pdb_parser.add_subparsers(dest="pdb_command", help="PDB commands")
    add_pdb_filter_subparser(pdb_subparsers)
    add_pdb_merge_subparser(pdb_subparsers)
    add_pdb_resname_subparser(pdb_subparsers)
    add_pdb_water_subparser(pdb_subparsers)
    add_pdb_unify_numbering_subparser(pdb_subparsers)

    # Prep subcommand and its nested commands
    prep_parser = subparsers.add_parser(
        "prep", help="Commands for preparing simulations"
    )
    parser._prep_parser = prep_parser  # type: ignore
    prep_subparsers = prep_parser.add_subparsers(
        dest="prep_command", help="PDB commands"
    )
    add_slurm_subparser(prep_subparsers)
    add_topo_subparser(prep_subparsers)

    return parser


def print_help_if_no_subcommand(
    parser: argparse.ArgumentParser, args: argparse.Namespace
) -> None:
    """Prints the help message of the appropriate parser if no subcommand is provided.

    This function checks if the user has provided a subcommand at the top level.
    If not, it prints the help message for the main parser and exits. It also
    handles the case where the user has selected the 'structure' top-level command
    but has not provided any nested subcommand under it, in which case it prints
    the help message for the 'structure' subparser. This ensures that users are
    provided with guidance on how to use the CLI when they don't provide a specific
    action to perform.

    Args:
        parser: The main argument parser object created
            by the `create_parser` function.
        args: The namespace object containing the parsed arguments.

    Raises:
        SystemExit: If no subcommand is provided, the function will call `sys.exit(0)`
            after printing the help message.

    Notes:
        The function checks the `args.command` attribute to determine if a top-level
        subcommand was provided. For the 'structure' command, it further checks for
        the presence of the `structure_command` attribute. The function relies on
        the `_structure_parser` attribute being attached to the main parser during
        its creation in the `create_parser` function to access the 'structure'
        subparser for printing its help message.
    """
    # Check for top-level subcommand
    if not args.command:
        parser.print_help()
        sys.exit(0)

    # Check for nested subcommand under 'structure'
    if args.command == "structure" and not getattr(args, "structure_command", None):
        structure_parser = getattr(parser, "_structure_parser", None)
        if structure_parser:
            structure_parser.print_help()
        else:
            parser.print_help()
        sys.exit(0)

    if args.command == "pdb" and not getattr(args, "pdb_command", None):
        pdb_parser = getattr(parser, "_pdb_parser", None)
        if pdb_parser:
            pdb_parser.print_help()
        else:
            parser.print_help()
        sys.exit(0)

    if args.command == "prep" and not getattr(args, "prep_command", None):
        prep_parser = getattr(parser, "_prep_parser", None)
        if prep_parser:
            prep_parser.print_help()
        else:
            parser.print_help()
        sys.exit(0)
