from simlify.utils import write_in_pdb_line


def write_atom_id(line: str, atom_id: int) -> str:
    r"""Unify residue ID in the PDB line based on previous ones.

    Args:
        line: Line that we are modifying.

    Returns:
        PDB line with the new atom ID.
    """
    line_start = 6
    line_stop = 12
    atom_line = str(atom_id).rjust(line_stop - line_start - 1) + " "
    new_line = write_in_pdb_line(line, atom_line, line_start, line_stop)
    return new_line
