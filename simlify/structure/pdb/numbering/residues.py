from loguru import logger

from simlify.utils import parse_resid, replace_in_pdb_line


def assign_resid(
    line: str, current_resid: str | None, prev_original_resid: str
) -> tuple[str, str]:
    r"""Determines residue ID based on a consistent numbering scheme.

    Args:
        line: Line that we are determining the residue ID to have.
        current_resid: Current residue ID that we are using.
        prev_original_resid: Previous residue ID from the PDB file that we are
            grouping together.

    Returns:
        Assigned residue ID for this line.

        Next `prev_original_resid`.
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
    r"""Unify residue ID in the PDB line based on previous ones.

    Args:
        line: Line that we are modifying.
        current_resid: Current residue ID that we are using.
        prev_original_resid: Original residue ID from previous line in the PDB file.

    Returns:
        PDB line with the new residue ID.

        Residue ID we are assigning this atom.

        The original residue ID from `line` (i.e., the next `prev_original_resid`).
    """
    assigned_resid, next_original_resid = assign_resid(
        line, current_resid, prev_original_resid
    )
    line_start = 22
    line_stop = 27
    new_line = replace_in_pdb_line(
        line,
        next_original_resid,
        assigned_resid.rjust(line_stop - line_start - 1) + " ",
        line_start,
        line_stop,
    )
    return new_line, assigned_resid, next_original_resid
