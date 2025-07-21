# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

"""
Module for assigning and unifying residue identifiers in PDB file lines to ensure
consistent numbering.
"""

from loguru import logger

from simlify.structure.pdb.utils import parse_resid, replace_in_pdb_line


def assign_resid(
    line: str, current_resid: str | None, prev_original_resid: str
) -> tuple[str, str]:
    r"""Assigns a consistent residue ID based on the sequence of residues encountered.

    This function takes a PDB line, the currently assigned residue ID, and the original
    residue ID from the previous line. It determines the residue ID for the current line
    based on whether the original residue ID has changed. If it's the first residue
    encountered (`current_resid` is None), it assigns the `prev_original_resid`. If the
    current line's original residue ID is the same as the previous one, it retains the
    `current_resid`. If the original residue ID has changed, it increments the
    `current_resid`.

    Args:
        line: The PDB line for which the residue ID is to be determined.
        current_resid: The residue ID that has been assigned to the previous
            residues. If this is the first residue, it will be `None`.
        prev_original_resid : The original residue ID as parsed from the previous
            PDB line. This is used to detect changes in the residue sequence.

    Returns:
        A tuple containing:

            -   The assigned residue ID for the current line, based on the consistent
                numbering scheme.
            -   The original residue ID parsed from the current `line`, which will serve
                as the `prev_original_resid` for the next line.

    Examples:
        >>> line1 = "ATOM      1  N   MET A   1      10.000  20.000  30.000  1.00 20.00           N"
        >>> line2 = "ATOM      2  CA  MET A   1      11.000  21.000  31.000  1.00 20.00           C"
        >>> line3 = "ATOM      3  C   ALA A   2      12.000  22.000  32.000  1.00 20.00           C"
        >>> assigned_id1, next_orig1 = assign_resid(line1, None, "   1")
        >>> print(f"{assigned_id1=}, {next_orig1=}")
        assigned_id1='   1', next_orig1='   1'
        >>> assigned_id2, next_orig2 = assign_resid(line2, assigned_id1, next_orig1)
        >>> print(f"{assigned_id2=}, {next_orig2=}")
        assigned_id2='   1', next_orig2='   1'
        >>> assigned_id3, next_orig3 = assign_resid(line3, assigned_id2, next_orig2)
        >>> print(f"{assigned_id3=}, {next_orig3=}")
        assigned_id3='2', next_orig3='   2'
    """
    next_original_resid = parse_resid(line).strip()
    logger.trace("Parsed residue ID from line: {}", next_original_resid)

    # We have our first residue.
    if current_resid is None:
        logger.debug("Current residue ID is None; must be our first atom.")
        assigned_resid = prev_original_resid
    else:
        # If the line's residue id is the same as the current original, then we should
        # group this atom with the previous one.
        assigned_resid = current_resid
        if next_original_resid != prev_original_resid:
            logger.trace("Parsed residue ID is not the same as previous.")
            logger.trace("Previous residue ID: {}", assigned_resid)
            assigned_resid = str(int(assigned_resid) + 1)
    logger.trace("Assigning residue ID: {}", assign_resid)
    return assigned_resid, next_original_resid


def unify_resid(
    line: str, current_resid: str | None, prev_original_resid: str
) -> tuple[str, str, str]:
    r"""Unifies the residue ID in a PDB line based on a consistent numbering scheme.

    This function takes a PDB line, the currently assigned residue ID, and the original
    residue ID from the previous line. It calls the `assign_resid` function to determine
    the unified residue ID for the current line and then replaces the original residue ID
    in the PDB line with this new, unified ID.

    Args:
        line: The PDB line to be modified.
        current_resid: The residue ID that has been assigned to the previous
            residues. If this is the first residue, it will be `None`.
        prev_original_resid: The original residue ID as parsed from the previous
            PDB line.

    Returns:
        A tuple containing:
            -   The modified PDB line with the unified residue ID.
            -   The residue ID that was assigned to this atom.
            -   The original residue ID parsed from the input `line`, which should
                be used as `prev_original_resid` for the next line.

    Examples:
        >>> line1 = "ATOM      1  N   MET A   1      10.000  20.000  30.000  1.00 20.00           N"
        >>> line2 = "ATOM      2  CA  MET A   1      11.000  21.000  31.000  1.00 20.00           C"
        >>> line3 = "ATOM      3  C   ALA A   2      12.000  22.000  32.000  1.00 20.00           C"
        >>> unified_line1, assigned_id1, next_orig1 = unify_resid(line1, None, "   1")
        >>> print(unified_line1.strip())
        ATOM      1  N   MET A   1      10.000  20.000  30.000  1.00 20.00           N
        >>> unified_line2, assigned_id2, next_orig2 = unify_resid(
        ...     line2, assigned_id1, next_orig1
        ... )
        >>> print(unified_line2.strip())
        ATOM      2  CA  MET A   1      11.000  21.000  31.000  1.00 20.00           C
        >>> unified_line3, assigned_id3, next_orig3 = unify_resid(
        ...     line3, assigned_id2, next_orig2
        ... )
        >>> print(unified_line3.strip())
        ATOM      3  C   ALA A   2      12.000  22.000  32.000  1.00 20.00           C
    """
    assigned_resid, next_original_resid = assign_resid(
        line, current_resid, prev_original_resid
    )
    line_start = 21
    line_stop = 26
    new_line = replace_in_pdb_line(
        line,
        next_original_resid,
        assigned_resid.rjust(line_stop - line_start),
        line_start,
        line_stop,
    )
    return new_line, assigned_resid, next_original_resid
