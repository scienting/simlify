from simlify.structure.pdb.utils import run_merge_pdbs


def add_pdb_merge_subparser(subparsers):
    r"""Command-line interface for merging PDB files"""
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


def cli_merge_pdbs(args, parser):
    r"""Command-line interface for merging PDB files"""
    if args.output is None:
        raise RuntimeError("--output must be specified")
    run_merge_pdbs(*args.pdb_paths, output_path=args.output)
