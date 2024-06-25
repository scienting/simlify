"""Standardizes residue ID numbering"""

import argparse
import os
from collections.abc import Iterable

from loguru import logger

from ..utils import parse_resid
from .atoms import write_atom_id
from .residues import unify_resid


def run_unify_numbering(
    pdb_path: str, output_path: str | None = None, reset_initial_resid: bool = True
) -> Iterable[str]:
    r"""Unify atom and residue numbering in the PDB lines.

    Args:
        pdb_path: Path to PDB file.
        output_path: Path to save new PDB file. If `None`, then no file is written.

    Returns:
        PDB file lines.
    """
    logger.info("Unify residue IDs from {}", os.path.abspath(pdb_path))
    with open(pdb_path, "r", encoding="utf-8") as f:
        pdb_lines: list[str] = f.readlines()

    # current_resid keeps track of our residue ID while parsing.
    # If this becomes None, then we restart our numbering at 1
    # (determined in assign_resid function).
    current_resid: str | None = None
    # current_chain is essentially the same thing as current_resid, but for the chain
    # ID.
    current_chain: str | None = None
    # parse_structure is used to turn parsing on or off depending on if the current
    # line contains any atomic information.
    parse_structure = False

    # atom_id keeps track of the atom serial number.
    atom_id: int = 1
    for i, line in enumerate(pdb_lines):
        logger.trace("Processing line number: {}", i)

        # When we hit a TER statement, we need to keep track/handle our IDs for the next
        # section.
        if line.startswith("TER"):
            # Catches and fixes duplicate TER statements
            if pdb_lines[i - 1].strip() == "TER":
                logger.debug("Fixing duplicate `TER` statements.")
                pdb_lines[i] = ""
                continue

            logger.debug("Encountered `TER`, storing ID information.")

            parse_structure = False
            logger.trace("Set parse_structure to False")

            current_resid = str(parse_resid(pdb_lines[i - 1]))
            logger.trace(f"Setting current_resid to {current_resid} (from last one)")

            if reset_initial_resid and current_chain is not None:
                current_chain = chr(ord(current_chain) + 1)
                logger.trace(f"Increasing current_chain to {current_chain}")

            pdb_lines[i] = "TER\n"
            continue

        # If we hit an ENDMDL statement, we are essentially providing another unique
        # PDB structure. So we reset everything.
        if line.startswith("ENDMDL"):
            logger.debug("Encountered `ENDMDL`; resetting all ID information")
            parse_structure = False
            current_resid = None
            current_chain = None
            atom_id = 1

        if line.startswith(("ATOM", "HETATM")):
            # Handle initializing the chain information
            chain_id = line[21]
            if current_chain is None:
                current_chain = chain_id

            # Activate coordinate parsing on first instance of ATOM or HETATM
            if not parse_structure:
                parse_structure = True

                # When we turn on parsing, this line contains the first atom
                # information. This is where we can reset our
                if not reset_initial_resid:
                    current_original_resid = str(parse_resid(line))
                else:
                    current_original_resid = "   1"

            # Write atom index.
            line = write_atom_id(line, atom_id)
            atom_id += 1

            # Write residue index.
            line, current_resid, current_original_resid = unify_resid(
                line, current_resid, current_original_resid
            )

            # Write chain.
            s_list = list(line)
            s_list[21] = current_chain
            line = "".join(s_list)

            # Update line
            pdb_lines[i] = line

    if output_path is not None:
        logger.info("Writing PDB file to {}", os.path.abspath(output_path))
        with open(output_path, "w+", encoding="utf-8") as f:
            f.writelines(pdb_lines)

    return pdb_lines


def cli_unify_numbering() -> None:
    r"""Command-line interface for unifying atom and residue numbering in PDB files."""
    parser = argparse.ArgumentParser(description="Unify residue IDs in PDB")
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
        "--keep_init_resid",
        action="store_false",
        help="Do not reset first residue ID to 1.",
    )

    args = parser.parse_args()
    run_unify_numbering(args.pdb_path, args.output, args.keep_init_resid)
