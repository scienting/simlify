# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

"""Module for writing or modifying atom identifiers within PDB file lines."""

from simlify.structure.pdb.utils import write_in_pdb_line


def write_atom_id(line: str, atom_id: int) -> str:
    r"""Writes a new atom ID into a specific PDB line.

    This function takes a PDB line and an integer representing the new atom ID.
    It formats the atom ID to fit within the standard atom serial number columns
    (columns 7-11, inclusive) of a PDB file line and then uses the `write_in_pdb_line`
    utility function to insert this new ID into the line.

    Args:
        line: The PDB line to be modified. This should be a standard PDB
            format line, typically starting with "ATOM" or "HETATM".
        atom_id: The new atom ID (serial number) to be written into the PDB line.

    Returns:
        The modified PDB line with the new atom ID written into the appropriate columns.

    Examples:
        >>> line = "ATOM      1  N   MET A   1      10.000  20.000  30.000  1.00 20.00           N"
        >>> new_line = write_atom_id(line, 100)
        >>> print(new_line)
        ATOM    100  N   MET A   1      10.000  20.000  30.000  1.00 20.00           N
    """
    line_start = 6  # zero-based index
    line_stop = 11
    field_width = line_stop - line_start  # 6 characters total
    # Format as right-justified in field_width, which includes the trailing space
    atom_line = str(atom_id).rjust(field_width)
    new_line = write_in_pdb_line(line, atom_line, line_start, line_stop)
    return new_line
