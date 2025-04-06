"""Standardizes residue ID numbering in PDB files."""

import argparse
import os
from collections.abc import Iterable

from loguru import logger

from simlify.structure.pdb.numbering.atoms import write_atom_id
from simlify.structure.pdb.numbering.residues import unify_resid
from simlify.structure.utils import parse_resid


def run_unify_numbering(
    pdb_path: str, output_path: str | None = None, reset_initial_resid: bool = True
) -> Iterable[str]:
    r"""Unifies the atom and residue numbering within a PDB file, ensuring sequential
    and consistent IDs.

    This function reads a PDB file and renumbers the atom serial numbers and residue IDs
    to be sequential, starting from 1 for the first atom and the first residue encountered.
    It also handles chain identifiers, incrementing them upon encountering a "TER" record
    if `reset_initial_resid` is True. Duplicate "TER" statements are removed, and "ENDMDL"
    records trigger a reset of atom and residue numbering, as well as the chain identifier.

    Args:
        pdb_path: The path to the input PDB file.
        output_path: The path to save the new PDB file with unified
            numbering. If `None`, no file is written, and the modified PDB lines are returned.
            Defaults to `None`.
        reset_initial_resid: If `True` (default), the residue numbering will
            start from 1 for the first residue in each chain. If `False`, the initial
            residue ID will be based on the original numbering in the PDB file for the
            first chain, and subsequent chains will continue sequentially.

    Returns:
        An iterable of strings, where each string is a line from the PDB file
            with the unified atom and residue numbering.

    Raises:
        FileNotFoundError: If the specified `pdb_path` does not exist.
        IOError: If there is an error reading the PDB file or writing to the output file.

    Notes:
        -   The function iterates through the PDB lines, tracking the current residue
            and chain IDs.
        -   When a "TER" record is encountered, it signifies the end of a chain, and
            the chain ID is incremented if `reset_initial_resid` is True.
        -   "ENDMDL" records indicate the start of a new model, and all numbering is
            reset.
        -   Atom serial numbers are simply incremented sequentially.
        -   Residue IDs are unified within each chain, potentially resetting to 1 at
            the start of a new chain.

    Examples:
        To unify the numbering in "input.pdb" and save it to "output.pdb":

        >>> unified_lines = run_unify_numbering("input.pdb", output_path="output.pdb")

        To unify the numbering but keep the initial residue ID of the first chain:

        >>> unified_lines = run_unify_numbering("input.pdb", output_path="output.pdb", reset_initial_resid=False)

        To unify the numbering and only get the lines without saving to a file:

        >>> unified_lines = run_unify_numbering("input.pdb")
        >>> for line in unified_lines:
        ...     print(line.strip())
        ...
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
    r"""Command-line interface for unifying atom and residue numbering in PDB files.

    This function provides a command-line tool to execute the `run_unify_numbering`
    functionality. It uses `argparse` to define the necessary command-line arguments
    for specifying the input PDB file, the output path for the modified PDB file,
    and an option to control whether the initial residue ID should be reset to 1.

    The command-line usage is as follows:

    ```bash
    python your_script_name.py input.pdb --output unified.pdb --keep_init_resid
    ```

    Where:

    -   `input.pdb` is the path to the PDB file to be processed.
    -   `--output unified.pdb` specifies the path for the output
        PDB file (optional).
    -   `--keep_init_resid` is an optional flag. If present, the initial residue
        ID will not be reset to 1.

    Raises:
        SystemExit: If the command-line arguments are invalid or if help is requested.

    See Also:
        `run_unify_numbering`: The underlying function that performs the PDB numbering
        unification.
    """
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
