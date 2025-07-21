# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

import argparse
import os
from collections.abc import Iterable

import MDAnalysis as mda
import numpy as np
from loguru import logger
from MDAnalysis import transformations as trans


def write_in_pdb_line(line: str, new: str, start: int | None, stop: int | None) -> str:
    r"""General function to write a new string into a specific portion of a PDB line.

    This function takes a PDB line and replaces a segment of it with a provided new string.
    It offers precise control over which part of the line is modified using start and stop indices.

    Args:
        line: The original PDB line to be modified.
        new: The new string to be inserted into the PDB line. This string should be
            formatted to match the expected width of the replaced segment, including
            any necessary spaces. For example, to represent the number 42 in a field
            that typically occupies 5 characters, the `new` string should be `"   42"`.
        start: The starting index (inclusive) of the slice in the `line` to be
            replaced. If `None`, the replacement starts from the beginning of the line.
        stop: The stopping index (exclusive) of the slice in the `line` to be
            replaced. If `None`, the replacement continues to the end of the line.

    Returns:
        The modified PDB line with the specified segment replaced by the `new` string.

    Examples:
        >>> line = "ATOM      1  N   MET A   1      10.000  20.000  30.000  1.00 20.00           N"
        >>> new_line = write_in_pdb_line(line, "   42", 6, 11)
        >>> print(new_line)
        ATOM     42  N   MET A   1      10.000  20.000  30.000  1.00 20.00           N
    """
    line = line[:start] + new + line[stop:]
    return line


def replace_in_pdb_line(
    line: str, orig: str, new: str, start: int | None, stop: int | None
) -> str:
    r"""General function to replace an original string with a new string within a
    specific portion of a PDB line.

    This function searches for a specific `orig` string within a defined slice of a
    PDB line. If the `orig` string is found, it is replaced with the provided `new`
    string. The replacement is constrained to the segment of the line specified by
    the `start` and `stop` indices.

    Args:
        line: The original PDB line to be examined and potentially modified.
        orig: The original string to search for within the specified slice of the
            `line`.
        new: The new string to replace the `orig` string if it is found. This string
            should be formatted to match the expected width of the replaced segment,
            including any necessary spaces. For example, to represent the residue
            number 42, the `new` string should be `"   42"` if the residue number
            field occupies 5 characters.
        start: The starting index (inclusive) of the slice in the `line` to be
            searched and where the replacement will occur. If `None`, the search
            and replacement start from the beginning of the line.
        stop: The stopping index (exclusive) of the slice in the `line` to be
            searched and where the replacement will occur. If `None`, the search
            and replacement continue to the end of the line.

    Returns:
        The modified PDB line where the `orig` string has been replaced by the
            `new` string within the specified slice, or the original line if `orig`
            was not found.

    Examples:
        >>> line = "ATOM      1  N   MET A   1      10.000  20.000  30.000  1.00 20.00           N"
        >>> new_line = replace_in_pdb_line(line, "MET", "ALA", 17, 20)
        >>> print(new_line)
        ATOM      1  N   ALA A   1      10.000  20.000  30.000  1.00 20.00           N
    """
    line_slice = line[start:stop]
    logger.trace("Slice gives us: '{}'", line_slice)
    if orig in line_slice:
        line_slice = new
    return line[:start] + line_slice + line[stop:]


def parse_resid(line: str) -> str:
    r"""Extracts the residue ID from a standard PDB file line.

    This function assumes the input `line` adheres to the standard PDB format for ATOM
    or HETATM records, where the residue ID is typically located in columns
    23-30 (inclusive).

    Args:
        line: A line from a PDB file that starts with either "ATOM" or "HETATM".

    Returns:
        The residue ID extracted from the line. This will typically include the residue
            sequence number and optionally an insertion code.

    Raises:
        IndexError: If the input `line` is shorter than 30 characters, accessing the
            residue ID slice will result in an `IndexError`.

    Examples:
        >>> line = "ATOM      1  N   MET A   1      10.000  20.000  30.000  1.00 20.00           N"
        >>> resid = parse_resid(line)
        >>> print(resid)
        '      1'
    """
    return line[22:30]


def parse_resname(line: str) -> str:
    r"""Extracts the residue name from a standard PDB file line.

    This function assumes the input `line` adheres to the standard PDB format for ATOM
    or HETATM records, where the residue name is typically located in columns
    18-21 (inclusive).

    Args:
        line: A line from a PDB file that starts with either "ATOM" or "HETATM".

    Returns:
        The residue name (e.g., "MET", "ALA", "HOH") extracted from the line.

    Raises:
        IndexError: If the input `line` is shorter than 21 characters, accessing the
            residue name slice will result in an `IndexError`.

    Examples:
        >>> line = "ATOM      1  N   MET A   1      10.000  20.000  30.000  1.00 20.00           N"
        >>> resname = parse_resname(line)
        >>> print(resname)
        'MET'
    """
    return line[17:21]


def parse_atomname(line: str) -> str:
    r"""Extracts the atom name from a standard PDB file line.

    This function assumes the input `line` adheres to the standard PDB format for ATOM
    or HETATM records, where the atom name is typically located in columns 14-17 (inclusive).

    Args:
        line: A line from a PDB file that starts with either "ATOM" or "HETATM".

    Returns:
        The atom name (e.g., "N", "CA", "O") extracted from the line.

    Raises:
        IndexError: If the input `line` is shorter than 17 characters, accessing the
            atom name slice will result in an `IndexError`.

    Examples:
        >>> line = "ATOM      1  N   MET A   1      10.000  20.000  30.000  1.00 20.00           N"
        >>> atomname = parse_atomname(line)
        >>> print(atomname)
        'N'
    """
    return line[13:17]


def keep_lines(
    lines: Iterable[str],
    record_types: tuple[str, ...] = ("ATOM", "HETATM", "TER", "END", "MODEL", "ENDMDL"),
) -> list[str]:
    r"""Filters a list of PDB file lines, retaining only those that start with specified
    record types.

    This function iterates through a given iterable of strings, which are assumed to be
    lines from a PDB file. It checks if each line begins with any of the record types
    provided in the `record_types` tuple. Only the lines that match one of these record
    types are included in the returned list.

    Args:
        lines: An iterable (e.g., a list) of strings, where each string represents
            a line from a PDB file.
        record_types: A tuple of strings representing the PDB record types to be kept.
            The default value is `("ATOM", "HETATM", "TER", "END", "MODEL", "ENDMDL")`,
            which includes the most common record types.

    Returns:
        A new list containing only the lines from the input `lines` that start with
            one of the specified `record_types`. The order of the lines in the output
            list will be the same as in the input.

    Examples:
        >>> pdb_lines = [
        ...     "HEADER    TITLE                                                   04-APR-25   NONE",
        ...     "ATOM      1  N   MET A   1      10.000  20.000  30.000  1.00 20.00           N",
        ...     "HELIX    1   1 MET A    1  THR A   4  1                                    4",
        ...     "HETATM  999  O   HOH     1      15.000  25.000  35.000  1.00 20.00           O",
        ...     "TER     1     MET A   1",
        ...     "END",
        ... ]
        >>> kept_lines = keep_lines(
        ...     pdb_lines, record_types=("ATOM", "HETATM", "TER", "END")
        ... )
        >>> for line in kept_lines:
        ...     print(line.strip())
        ATOM      1  N   MET A   1      10.000  20.000  30.000  1.00 20.00           N
        HETATM  999  O   HOH     1      15.000  25.000  35.000  1.00 20.00           O
        TER     1     MET A   1
        END
    """
    logger.info("Keeping the following record types: {}", ", ".join(record_types))
    return [line for line in lines if line.startswith(record_types)]


def run_filter_pdb(
    pdb_path: str,
    output_path: str | None = None,
    record_types: tuple[str, ...] | None = None,
) -> list[str]:
    r"""Reads a PDB file and keeps only the lines that contain specified record types.

    This function takes the path to a PDB file, reads its contents, filters the lines
    to retain only those that start with the record types specified in the
    `record_types` argument, and optionally writes the filtered lines to a new PDB file.

    Args:
        pdb_path: The path to the input PDB file.
        output_path: The path to the output PDB file where the filtered lines will be
            written. If `None`, no new file is created, and the filtered lines are
            only returned.
        record_types: A tuple of strings representing the PDB record types to be kept.
            If `None`, the default record types
            `("ATOM", "HETATM", "TER", "END", "MODEL", "ENDMDL")` are used.

    Returns:
        A list of strings, where each string is a line from the PDB file that starts
            with one of the specified `record_types`.

    Raises:
        FileNotFoundError: If the `pdb_path` does not exist.
        IOError: If there is an error reading or writing the PDB files.

    Examples:
        To filter a PDB file named "input.pdb" and save the result to "output.pdb",
        keeping only ATOM and TER records:

        >>> run_filter_pdb("input.pdb", "output.pdb", record_types=("ATOM", "TER"))

        To filter a PDB file and only get the lines without writing to a new file:

        >>> filtered_lines = run_filter_pdb(
        ...     "input.pdb", record_types=("ATOM", "HETATM")
        ... )
        >>> for line in filtered_lines:
        ...     print(line.strip())
    """
    logger.info("Filtering PDB lines of {}", os.path.abspath(pdb_path))
    with open(pdb_path, "r", encoding="utf-8") as f:
        pdb_lines: list[str] = f.readlines()

    if record_types is None:
        record_types = ("ATOM", "HETATM", "TER", "END", "MODEL", "ENDMDL")
    out_lines = keep_lines(pdb_lines, record_types)

    if output_path is not None:
        logger.info("Writing PDB file to {}", os.path.abspath(output_path))
        with open(output_path, "w+", encoding="utf-8") as f:
            f.writelines(out_lines)

    return out_lines


def run_merge_pdbs(*pdb_paths: str, output_path: str | None = None) -> mda.Universe:
    r"""Merges multiple PDB files into a single MDAnalysis Universe object.

    This function takes a variable number of PDB file paths as input. It loads the
    first PDB file into an MDAnalysis Universe object and then iteratively adds the
    atoms from the subsequent PDB files. It assumes that the residue indices are
    consistent across all input PDB files. The merging process attempts to add missing
    atom types to existing residues based on the information in the later PDB files.
    Duplicate atoms (based on their coordinates) are removed, and the atoms within each
    residue are sorted by their type. Finally, some topology attributes
    that might interfere with other programs are removed.

    Args:
        *pdb_paths: A variable number of strings, where each string is the path to a PDB
            file. The order of the paths is important, as the first PDB file sets the
            initial structure, and subsequent files are used to add missing atoms.
        output_path: The path to save the merged PDB structure to a new
            file. If `None`, no file is written.

    Returns:
        An MDAnalysis Universe object containing all the atoms from the input PDB
            files, with duplicate atoms removed and atoms within each residue sorted.

    Raises:
        FileNotFoundError: If any of the provided `pdb_paths` do not exist.
        IOError: If there is an error reading any of the PDB files.

    Notes:
        -   This function prioritizes the atom information from the PDB files provided
            later in the argument list when resolving missing atom types within a
            residue.
        -   The function removes the "segids" topology attribute from the merged
            Universe, as it can sometimes cause issues with programs like pdb4amber.

    Examples:
        To merge two PDB files, "file1.pdb" and "file2.pdb", and save the result to
        "merged.pdb":

        >>> merged_universe = run_merge_pdbs(
        ...     "file1.pdb", "file2.pdb", output_path="merged.pdb"
        ... )

        To merge multiple PDB files without saving to a new file:

        >>> merged_universe = run_merge_pdbs("file1.pdb", "file2.pdb", "file3.pdb")
    """
    logger.info("Merging PDB files")
    u = mda.Universe(pdb_paths[0])
    atoms_to_add = [mda.Universe(pdb_path).atoms for pdb_path in pdb_paths[1:]]
    u_to_add = mda.core.universe.Merge(*atoms_to_add)

    for residue in u_to_add.residues:
        logger.debug("Processing residue ID of {}", residue.resid)
        available_types = {atom.name for atom in residue.atoms}
        logger.trace("Available atom types: {}", available_types)
        present_types = {atom.name for atom in u.select_atoms(f"resid {residue.resid}")}
        logger.trace("Current atom types: {}", present_types)
        missing_types = {
            atype for atype in available_types if atype not in present_types
        }
        if len(missing_types) > 0:
            logger.trace("Missing atom types: {}", missing_types)
            add_atoms = [atom for atom in residue.atoms if atom.name in missing_types]
            u = mda.core.universe.Merge(u.atoms, mda.AtomGroup(add_atoms))

    # Get the indices of duplicate atoms
    coordinates = u.atoms.positions
    _, unique_indices = np.unique(coordinates, axis=0, return_index=True)
    if len(unique_indices) < len(coordinates):
        logger.info("Cleaning up duplicate atoms")
        duplicate_indices = np.setdiff1d(np.arange(len(coordinates)), unique_indices)
        u.atoms = u.atoms[
            np.isin(np.arange(len(coordinates)), duplicate_indices, invert=True)
        ]

    # Group atoms by residue
    residue_groups = u.atoms.groupby("resids")
    u = mda.core.universe.Merge(
        *[atoms.sort("types") for _, atoms in residue_groups.items()]
    )

    # Remove some Merge artifacts that messes with pdb4amber
    u.del_TopologyAttr("segids")

    if output_path is not None:
        logger.info("Writing merged PDB at {}", output_path)
        u.atoms.write(output_path)
    return u.atoms


def run_write_pdb(
    file_paths: Iterable[str],
    output_path: str,
    selection_str: str | None = None,
    stride: int = 1,
) -> None:
    r"""Writes a PDB file from a set of topology and coordinate files,
    potentially applying a selection and stride.

    This function takes a list of file paths that can be read by MDAnalysis to create
    a Universe object. It then iterates through the trajectory of this Universe,
    and at each time step (optionally with a specified stride), it writes the
    coordinates of the selected atoms to a PDB file.

    Args:
        file_paths: An iterable of strings, where each string is the path to a
            topology or coordinate file that can be loaded by MDAnalysis
            (e.g., a topology file like a PRMTOP or a coordinate file like a TRR, DCD,
            or PDB). If multiple files are provided, the first is typically the
            topology, and the rest are coordinate files.
        output_path: The path to the output PDB file that will be created or
            overwritten.
        selection_str: An MDAnalysis selection string that specifies which atoms to
            write to the PDB file. If `None`, all atoms in the current frame are
            written.
        stride: An integer specifying the stride for writing frames from the
            trajectory. Only frames where the frame number modulo `stride` is 0 will
            be written. A stride of 1 means every frame is written. Defaults to 1.

    Raises:
        FileNotFoundError: If any of the files specified in `file_paths` do not exist.
        IOError: If there is an error reading the input files or writing the output PDB file.
        ValueError: If the `file_paths` iterable is empty.

    Examples:
        To write all atoms from a TRR trajectory file "traj.trr" and topology file
        "top.pdb" to a PDB file "output.pdb":

        >>> run_write_pdb(["top.pdb", "traj.trr"], "output.pdb")

        To write only the protein atoms with a stride of 10:

        >>> run_write_pdb(
        ...     ["top.pdb", "traj.dcd"],
        ...     "protein.pdb",
        ...     selection_str="protein",
        ...     stride=10,
        ... )

        To write all atoms from a single PDB file to another PDB file (effectively copying it):

        >>> run_write_pdb(["input.pdb"], "output.pdb")
    """
    u = mda.Universe(*file_paths)
    with mda.Writer(output_path, multiframe=True) as W:
        for ts in u.trajectory:
            if ts.frame % stride == 0:
                if isinstance(selection_str, str):
                    atoms = u.select_atoms(selection_str)
                else:
                    atoms = u.atoms
                W.write(atoms)


def cli_write_pdb() -> None:
    r"""Command-line interface for writing a PDB file from topology and coordinate
    files.

    This function sets up an argument parser to allow users to write PDB files from
    the command line. It takes arguments for the output path, input files, an optional
    MDAnalysis selection string, and an optional stride for writing frames.

    The command-line usage is as follows:

    ```bash
    python your_script_name.py output.pdb --files top.pdb traj.dcd --select "protein and name CA" --stride 10
    ```

    This would write a PDB file named "output.pdb" containing only the alpha carbon
    atoms of the protein from the trajectory "traj.dcd" (with topology in "top.pdb"),
    taking every 10th frame.

    Raises:
        SystemExit: If the command-line arguments are invalid or if help is requested.

    See Also:
        `run_write_pdb`: The underlying function that performs the PDB writing.
    """
    parser = argparse.ArgumentParser(
        description="Write PDB from topology and coordinate files."
    )
    parser.add_argument(
        "output_path",
        type=str,
        nargs="?",
        help="PDB file to write",
    )
    parser.add_argument(
        "--files",
        type=str,
        nargs="+",
        help="Files to load into MDAnalysis universe.",
    )
    parser.add_argument(
        "--select",
        type=str,
        nargs="*",
        help="Selection string for MDAnalysis universe.",
    )
    parser.add_argument(
        "--stride",
        type=int,
        nargs="?",
        help="Stride of trajectory when writing.",
        default=1,
    )
    args = parser.parse_args()
    if args.select is not None:
        args.select = " ".join(args.select)
    run_write_pdb(args.files, args.output_path, args.select, args.stride)


def run_align_pdb(
    pdb_path: str,
    out_path: str,
    selection_str: str | None = None,
) -> None:
    r"""Aligns the structure within a PDB file to a reference configuration based on a
    selection of atoms.

    This function loads a PDB file into an MDAnalysis Universe, selects a subset of
    atoms based on the `selection_str`, and then performs a rigid-body fit of these
    atoms to their initial positions in the first frame of the trajectory.
    The transformation (rotation and translation) that achieves this fit
    is then applied to all atoms in all frames of the trajectory. The aligned
    trajectory is then written to a new PDB file.

    Args:
        pdb_path: The path to the input PDB file containing the structure to be aligned.
        out_path: The path to the output PDB file where the aligned structure will be
            written.
        selection_str: An MDAnalysis selection string that specifies the group
            of atoms to be used for the alignment. If `None`, all atoms in the
            structure are used for alignment.

    Raises:
        FileNotFoundError: If the input `pdb_path` does not exist.
        IOError: If there is an error reading the input PDB file or writing the output
            PDB file.
        ValueError: If the `selection_str` does not select any atoms.

    Notes:
        -   The alignment is performed against the conformation in the first frame of
            the input PDB file.
        -   This function is useful for removing overall translation and rotation from
            a structural ensemble.

    Examples:
        To align a PDB file "input.pdb" to its first frame using all atoms and save the
        result to "aligned.pdb":

        >>> run_align_pdb("input.pdb", "aligned.pdb")

        To align only the backbone atoms (N, CA, C) of the protein:

        >>> run_align_pdb(
        ...     "input.pdb",
        ...     "aligned_backbone.pdb",
        ...     selection_str="protein and backbone",
        ... )
    """
    u = mda.Universe(pdb_path)
    ag = u.select_atoms(selection_str)
    u_ref = u.copy()
    ag_ref = u_ref.select_atoms(selection_str)
    workflow = (trans.fit_rot_trans(ag, ag_ref),)
    u.trajectory.add_transformations(*workflow)
    with mda.Writer(out_path, multiframe=True) as W:
        for ts in u.trajectory:
            W.write(u.atoms)


def cli_align_pdb() -> None:
    r"""Command-line interface for aligning a PDB file.

    This function sets up an argument parser to allow users to align PDB files from
    the command line. It takes arguments for the input PDB path, the output PDB path,
    and an optional MDAnalysis selection string to specify which atoms should be used
    for the alignment.

    The command-line usage is as follows:

    ```bash
    python your_script_name.py input.pdb aligned.pdb --selection "protein and name CA"
    ```

    This would align the PDB file "input.pdb" to its first frame based on the alpha
    carbon atoms of the protein and save the aligned structure to "aligned.pdb".

    Raises:
        SystemExit: If the command-line arguments are invalid or if help is requested.

    See Also:
        `run_align_pdb`: The underlying function that performs the PDB alignment.
    """
    parser = argparse.ArgumentParser(description="Align PDB file to some selection.")
    parser.add_argument(
        "pdb_path",
        type=str,
        nargs="?",
        help="PDB file to load",
    )
    parser.add_argument(
        "out_path",
        type=str,
        nargs="?",
        help="PDB file to write",
    )
    parser.add_argument(
        "--selection",
        type=str,
        nargs="*",
        help="Selection string for MDAnalysis universe.",
    )
    args = parser.parse_args()
    if args.selection is not None:
        args.selection = " ".join(args.selection)
    run_align_pdb(args.pdb_path, args.out_path, args.selection)
