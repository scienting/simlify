"""
Structure module command-line interface.
"""

import argparse

from simlify.structure import extract_atoms


def cli_extract_atoms() -> None:
    r"""Command-line interface for extracting atoms or frames in structures."""
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
    extract_atoms(select, args.topo, args.output, args.coords)
