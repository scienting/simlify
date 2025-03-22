import argparse
import sys

from simlify.cli.structure.extract import add_extract_subparser


def create_parser() -> argparse.ArgumentParser:
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
    structure_subparsers = structure_parser.add_subparsers(
        dest="structure_command", help="Structure commands"
    )
    add_extract_subparser(structure_subparsers)

    return parser


def print_help_if_no_subcommand(
    parser: argparse.ArgumentParser, args: argparse.Namespace
) -> None:
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
