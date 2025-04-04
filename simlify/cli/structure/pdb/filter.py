from simlify.structure.pdb.utils import run_filter_pdb


def add_pdb_filter_subparser(subparsers):
    r"""Command-line interface for filtering structure."""
    parser = subparsers.add_parser("filter", description="Filter PDB lines")
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
        "--records",
        type=str,
        nargs="*",
        help="Records to keep in the PDB file.",
    )
    parser.set_defaults(func=lambda args: cli_filter_pdb(args, parser))
    return parser


def cli_filter_pdb(args, parser):
    r"""Command-line interface for filtering PDB file lines"""
    run_filter_pdb(args.pdb_path, args.output, args.records)
