"""Module for modifying lines within PDB files, including functionalities for replacing
atom and residue names.

This module provides a set of utility functions to manipulate the content of Protein
Data Bank (PDB) files. It includes functions to perform general line modifications
based on filtering criteria, as well as specific functions for replacing atom names
and residue names within the ATOM and HETATM records of a PDB file. Additionally, it
offers a function to standardize the atom names of water molecules.
"""

from typing import Any

import os
from collections.abc import Callable, Iterable

from loguru import logger

from simlify.structure.pdb.utils import (
    parse_atomname,
    parse_resid,
    parse_resname,
    replace_in_pdb_line,
)


def modify_lines(
    pdb_lines: Iterable[str],
    fn_process: Callable[[str, str, str, int | None, int], str],
    fn_args: Iterable[Any],
    fn_filter: Callable[[str], str] | None = None,
    include: list[str] | None = None,
    exclude: list[str] | None = None,
) -> list[str]:
    r"""General function to modify specific lines in a PDB file based on filtering.

    This function iterates through a list of PDB lines and applies a processing function
    (`fn_process`) to lines that meet certain criteria defined by an optional filter
    function (`fn_filter`) and inclusion/exclusion lists.

    Args:
        pdb_lines: An iterable of strings, where each string represents a line
            from a PDB file.
        fn_process: A callable function that takes a PDB line as its first argument,
            followed by the elements of `fn_args`, and returns a modified PDB line.
            This function is responsible for the actual modification of the line.
        fn_args: An iterable containing additional arguments to be passed to the
            `fn_process` function after the PDB line itself.
        fn_filter: An optional callable function that takes
            a PDB line as input and returns a string. This string is then used to
            check against the `include` and `exclude` lists. If `None`, all ATOM and
            HETATM lines are processed.
        include: An optional list of strings. If `fn_filter` is provided,
            only lines for which the result of `fn_filter` is present in this list
            will be processed by `fn_process`.
        exclude: An optional list of strings. If `fn_filter` is provided,
            lines for which the result of `fn_filter` is present in this list
            will *not* be processed by `fn_process`. Defaults to `None`.

    Returns:
        A list of modified PDB lines. Lines that did not meet the filtering criteria
            or were not ATOM or HETATM records are returned unchanged.

    Notes:
        - The `fn_filter` function should be designed to extract a specific piece of information
          from the PDB line (e.g., residue name, atom name) that can be used for inclusion or
          exclusion.
        - If both `include` and `exclude` are provided and a filtered value is present in both,
          the line will be processed if it's in `include`. Exclusion takes precedence if only
          `exclude` is provided.

    Examples:
        To replace "CA" atom names with "CB" only in residues named "GLY":

        >>> pdb_lines = [
        ...     "ATOM      1  CA  GLY A   1       ...",
        ...     "ATOM      2  CB  ALA A   2       ...",
        ... ]
        >>> def get_resname(line):
        ...     return parse_resname(line).strip()
        >>> modified = modify_lines(
        ...     pdb_lines,
        ...     replace_in_pdb_line,
        ...     ("CA ", "CB ", 13, 17),
        ...     fn_filter=get_resname,
        ...     include=["GLY"],
        ... )
        >>> for line in modified:
        ...     print(line)
        ATOM      1  CB  GLY A   1       ...
        ATOM      2  CB  ALA A   2       ...
    """
    modified_lines = []
    if callable(fn_filter):
        logger.debug("Filter function is provided.")
    for line in pdb_lines:
        if "ATOM" in line or "HETATM" in line:
            if callable(fn_filter):
                filter_return = fn_filter(line).strip()
                logger.trace("Filter function provided: {}", filter_return)
                # Line is in included and should be processed.
                if (include is not None) and (filter_return in include):
                    line = fn_process(line, *fn_args)
                # Line is not in excluded and should be processed.
                if (exclude is not None) and (filter_return not in exclude):
                    line = fn_process(line, *fn_args)
            else:
                line = fn_process(line, *fn_args)
            modified_lines.append(line)

    return modified_lines


def replace_atom_names(
    pdb_lines: Iterable[str], orig_atom_name: str, new_atom_name: str
) -> list[str]:
    r"""Replaces all occurrences of a specified original atom name with a new atom
    name in a list of PDB lines.

    This function iterates through the provided PDB lines and, for each ATOM or
    HETATM record, it checks if the atom name matches the `orig_atom_name`.
    If it does, the atom name is replaced with the `new_atom_name`. The atom names are
    stripped of leading/trailing whitespace and left-justified
    to a length of 4 characters to ensure proper formatting in the PDB file.

    Args:
        pdb_lines: An iterable of strings, where each string represents a line
            from a PDB file.
        orig_atom_name: The original atom name to be replaced.
        new_atom_name: The new atom name to replace the original one.

    Returns:
        list[str]: A list of PDB lines with the specified atom names replaced.

    Examples:
        >>> pdb_lines = [
        ...     "ATOM      1  CA  ALA A   1       ...",
        ...     "ATOM      2  CB  ALA A   1       ...",
        ... ]
        >>> modified_lines = replace_atom_names(pdb_lines, "CA", "CB")
        >>> for line in modified_lines:
        ...     print(line)
        ATOM      1  CB  ALA A   1       ...
        ATOM      2  CB  ALA A   1       ...
    """
    orig_atom_name = orig_atom_name.strip().ljust(4)
    new_atom_name = new_atom_name.strip().ljust(4)
    return modify_lines(
        pdb_lines, replace_in_pdb_line, (orig_atom_name, new_atom_name, 13, 17)
    )


def replace_residue_names(
    pdb_lines: Iterable[str],
    orig_resname: str,
    new_resname: str,
    fn_filter: Callable[[str], str] | None = None,
    include: list[str] | None = None,
    exclude: list[str] | None = None,
) -> list[str]:
    r"""Replaces all occurrences of a specified original residue name with a new
    residue name in a list of PDB lines.

    This function iterates through the provided PDB lines and, for each ATOM or
    HETATM record, it checks if the residue name matches the `orig_resname`.
    If it does, the residue name is replaced with the `new_resname`. The residue names
    are stripped of leading/trailing whitespace and left-justified to a length of 4
    characters to ensure proper formatting in the PDB file. Optionally, a filter
    function and inclusion/exclusion lists can be used to control which lines
    are processed.

    Args:
        pdb_lines: An iterable of strings, where each string represents a line
            from a PDB file.
        orig_resname: The original residue name to be replaced.
        new_resname: The new residue name to replace the original one.
        fn_filter: An optional callable function that takes a PDB line as input and
            returns a string (e.g., residue ID) for filtering.
        include: An optional list of strings. Only lines where the result of
            `fn_filter` is in this list will have their residue names replaced.
        exclude: An optional list of strings. Lines where the result of `fn_filter`
            is in this list will *not* have their residue names replaced.

    Returns:
        A list of PDB lines with the specified residue names replaced, subject to any
            provided filtering.

    Examples:
        To replace all "MET" residues with "ALA":

        >>> pdb_lines = [
        ...     "ATOM      1  N   MET A   1       ...",
        ...     "ATOM      2  CA  MET A   1       ...",
        ... ]
        >>> modified_lines = replace_residue_names(pdb_lines, "MET", "ALA")
        >>> for line in modified_lines:
        ...     print(line)
        ATOM      1  N   ALA A   1       ...
        ATOM      2  CA  ALA A   1       ...

        To replace "MET" with "ALA" only in residue ID "1":

        >>> pdb_lines = [
        ...     "ATOM      1  N   MET A   1       ...",
        ...     "ATOM      2  CA  MET A   1       ...",
        ...     "ATOM      3  C   MET A   2       ...",
        ... ]
        >>> def get_resid(line):
        ...     return parse_resid(line).strip()
        >>> modified_lines = replace_residue_names(
        ...     pdb_lines, "MET", "ALA", fn_filter=get_resid, include=["1"]
        ... )
        >>> for line in modified_lines:
        ...     print(line)
        ATOM      1  N   ALA A   1       ...
        ATOM      2  CA  ALA A   1       ...
        ATOM      3  C   MET A   2       ...
    """
    orig_resname = orig_resname.strip().ljust(4)
    new_resname = new_resname.strip().ljust(4)
    logger.info("Renaming '{}' to '{}'", orig_resname, new_resname)
    return modify_lines(
        pdb_lines,
        replace_in_pdb_line,
        (orig_resname, new_resname, 17, 21),
        fn_filter,
        include,
        exclude,
    )


def run_replace_resnames(
    pdb_path: str,
    resname_map: dict[str, str],
    output_path: str | None = None,
    fn_filter: Callable[[str], str] | None = None,
    include: list[str] | None = None,
    exclude: list[str] | None = None,
) -> list[str]:
    r"""Replaces multiple residue names in a PDB file based on a provided mapping.

    This function reads a PDB file, iterates through a dictionary that maps original
    residue names to new residue names, and applies the replacement using the
    `replace_residue_names` function for each mapping. The modified PDB lines are
    then either returned or written to a new file. Optional filtering based on a
    function and inclusion/exclusion lists can be applied during the replacement
    process for each residue name in the map.

    Args:
        pdb_path: The path to the input PDB file.
        resname_map: A dictionary where the keys are the original residue names
            to be replaced, and the values are the corresponding new residue names.
        output_path: The path to save the new PDB file with the replaced residue names.
            If `None`, no file is written, and the modified PDB lines are returned.
        fn_filter: An optional callable function that takes a PDB line as input and
            returns a string for filtering during each residue name replacement.
        include: An optional list of strings. Only lines where the result of
            `fn_filter` is in this list will have their residue names replaced for each
            mapping in `resname_map`. Defaults to `None`.
        exclude: An optional list of strings. Lines where the result of `fn_filter`
            is in this list will *not* have their residue names replaced for each
            mapping in `resname_map`.

    Returns:
        A list of PDB lines with the residue names replaced according to the
            `resname_map`, subject to any provided filtering.

    Raises:
        FileNotFoundError: If the specified `pdb_path` does not exist.
        IOError: If there is an error reading the PDB file or writing to the output
            file.

    Examples:
        To replace all "MET" residues with "ALA" and all "GLU" residues with "ASP" in
        "input.pdb" and save the result to "output.pdb":

        >>> resname_mapping = {"MET": "ALA", "GLU": "ASP"}
        >>> modified_lines = run_replace_resnames(
        ...     "input.pdb", resname_mapping, output_path="output.pdb"
        ... )

        To perform the same replacement but only for residues with ID "1":

        >>> def get_resid(line):
        ...     return parse_resid(line).strip()
        >>> resname_mapping = {"MET": "ALA", "GLU": "ASP"}
        >>> modified_lines = run_replace_resnames(
        ...     "input.pdb",
        ...     resname_mapping,
        ...     output_path="filtered_output.pdb",
        ...     fn_filter=get_resid,
        ...     include=["1"],
        ... )
    """
    logger.info("Renaming residue names {}", os.path.abspath(pdb_path))
    with open(pdb_path, "r", encoding="utf-8") as f:
        pdb_lines: list[str] = f.readlines()

    for orig_resname, new_resname in resname_map.items():
        pdb_lines = replace_residue_names(
            pdb_lines, orig_resname, new_resname, fn_filter, include, exclude
        )

    if output_path is not None:
        logger.info("Writing PDB file to {}", os.path.abspath(output_path))
        with open(output_path, "w+", encoding="utf-8") as f:
            f.writelines(pdb_lines)

    return pdb_lines


def run_unify_water_labels(
    pdb_path: str,
    atom_map: dict[str, str] | None = None,
    water_resname: str = "WAT",
    water_atomnames: dict[str, Iterable[str]] | None = None,
    output_path: str | None = None,
) -> Iterable[str]:
    r"""Ensures that water molecule atom names are consistently labeled as 'O', 'H1',
    and 'H2'.

    This function processes a PDB file to standardize the atom names of water
    molecules. It identifies water residues based on the `water_resname` and then
    renames their atoms to 'O' for oxygen and 'H1' and 'H2' for the two hydrogen atoms.
    The hydrogen atoms are assigned 'H1' and 'H2' based on their sequential appearance
    within each water residue in the PDB file.

    Args:
        pdb_path: The path to the input PDB file.
        atom_map: A dictionary mapping the standard water atom names ('O', 'H1', 'H2')
            to the desired names. If `None`, it defaults to
            `{'O': 'O', 'H1': 'H1', 'H2': 'H2'}`. This allows for customization of the
            output atom names if needed. Defaults to `None`.
        water_resname: The residue name used to identify water molecules in the PDB
            file.
        water_atomnames: A dictionary specifying the original atom names that should
            be considered as oxygen and hydrogen atoms of water. The keys should be
            'O' and 'H', and the values should be iterables of possible atom names.
            If `None`, it defaults to `{'O': ['OW'], 'H': ['HW']}`.
        output_path: The path to save the new PDB file with the unified water atom
            labels. If `None`, no file is written, and the modified PDB lines are
            returned.

    Returns:
        An iterable of PDB lines with the unified water atom labels.

    Raises:
        FileNotFoundError: If the specified `pdb_path` does not exist.
        IOError: If there is an error reading the PDB file or writing to the output
            file.

    Warning:
        This function has not been thoroughly tested and might not handle all edge
        cases correctly. Use with caution.

    Examples:
        To unify water atom labels in "input.pdb" using the default settings and save
        to "unified_water.pdb":

        >>> modified_lines = run_unify_water_labels(
        ...     "input.pdb", output_path="unified_water.pdb"
        ... )

        To specify a different water residue name and atom name mapping:

        >>> atom_mapping = {"O": "OXT", "H1": "HT1", "H2": "HT2"}
        >>> original_water_names = {"O": ["SOL"], "H": ["HY"]}
        >>> modified_lines = run_unify_water_labels(
        ...     "input.pdb",
        ...     atom_map=atom_mapping,
        ...     water_resname="SOL",
        ...     water_atomnames=original_water_names,
        ...     output_path="custom_water.pdb",
        ... )
    """
    logger.info("Renaming water atom names in {}", os.path.abspath(pdb_path))

    logger.debug("Water residue name: {}", water_resname)
    if atom_map is None:
        atom_map = {"O": "O", "H1": "H1", "H2": "H2"}
    logger.debug("O atom name: {}", atom_map["O"])
    logger.debug("H1 atom name: {}", atom_map["H1"])
    logger.debug("H2 atom name: {}", atom_map["H2"])
    if water_atomnames is None:
        water_atomnames = {"O": ["OW"], "H": ["HW"]}

    with open(pdb_path, "r", encoding="utf-8") as f:
        pdb_lines: list[str] = f.readlines()

    water_h_resids = []

    with open(pdb_path, "r", encoding="utf-8") as f:
        pdb_lines = f.readlines()

    for i, line in enumerate(pdb_lines):
        logger.trace("Working on: {}", line.strip())
        if parse_resname(line).strip() == water_resname:
            logger.trace("Line is water residue")
            original_atom_name = parse_atomname(line).strip()

            if original_atom_name in water_atomnames["O"]:
                logger.trace(
                    "Atom name, {}, matches an oxygen type", original_atom_name
                )
                pdb_lines[i] = replace_atom_names(
                    [line], original_atom_name, atom_map["O"]
                )[0]
            elif original_atom_name.startswith("H"):
                # Use the original hydrogen atom name as the key in the dictionary
                resid = parse_resid(line).strip()
                if resid not in water_h_resids:
                    water_h_resids.append(resid)
                    atom_map_name = atom_map["H1"]
                else:
                    atom_map_name = atom_map["H2"]
                pdb_lines[i] = replace_atom_names(
                    [line],
                    original_atom_name,
                    atom_map_name,
                )[0]

    if output_path is not None:
        logger.info("Writing PDB file to {}", os.path.abspath(output_path))
        with open(output_path, "w+", encoding="utf-8") as f:
            f.writelines(pdb_lines)

    return pdb_lines
