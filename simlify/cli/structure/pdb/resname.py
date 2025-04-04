from simlify.structure.pdb.names import run_replace_resnames
from simlify.structure.pdb.utils import parse_resid


def add_pdb_resname_subparser(subparsers):
    r"""Command-line interface for renaming residues."""
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


def cli_replace_resnames(args, parser):
    r"""Command-line interface for renaming residues."""
    resname_map = {args.current_resname: args.new_resname}
    if (args.include is not None) or (args.exclude is not None):
        fn_filter = parse_resid
    else:
        fn_filter = None
    run_replace_resnames(
        args.pdb_path, resname_map, args.output, fn_filter, args.include, args.exclude
    )
