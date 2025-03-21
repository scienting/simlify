"""
Structure module command-line interface.
"""

import argparse

from simlify.structure.select import run_select_atoms


def cli_select_atoms() -> None:
    r"""Command-line interface for selecting atoms in structures."""
    parser = argparse.ArgumentParser(description="Select atoms from PDB file")
    parser.add_argument(
        "topo",
        type=str,
        nargs="?",
        help="Path to topology file.",
    )
    parser.add_argument(
        "output",
        type=str,
        nargs="?",
        help="Path to coordinate file to save.",
    )
    parser.add_argument(
        "--select",
        type=str,
        nargs="+",
        help="MDAnalysis selection string.",
    )
    parser.add_argument(
        "--coords",
        type=str,
        nargs="*",
        help="Paths to coordinate files. Will concatenate these into one file.",
    )
    args = parser.parse_args()
    select = " ".join(args.select)
    run_select_atoms(select, args.topo, args.output, args.coords)
