import argparse

import MDAnalysis as mda
from MDAnalysis import transformations as trans


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
