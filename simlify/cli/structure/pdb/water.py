from simlify.structure.pdb.names import run_unify_water_labels


def add_pdb_water_subparser(subparsers):
    r"""Command-line interface for filtering structure."""
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


def cli_unify_water_labels(args, parser):
    r"""Command-line interface for unifying water residue and atom names."""
    args = parser.parse_args()
    run_unify_water_labels(args.pdb_path, output_path=args.output)
