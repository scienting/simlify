"""
Command-line interface for minimizing the simulation box size of a PDB file.
"""

from simlify.structure.pdb.orientation import run_minimize_box


def add_min_box_subparser(subparsers):
    r"""Adds the 'min-box' subcommand to the Simlify CLI for minimizing the box size.

    This function configures an `argparse` subparser named 'min-box' which allows
    users to specify a PDB file and an optional output path. The 'min-box'
    command will attempt to rotate the system within the PDB file to achieve a
    smaller bounding box volume.

    Args:
        subparsers: An `argparse._SubParsersAction` object where the 'min-box'
            subparser will be added. This is typically obtained by calling
            `add_subparsers()` on an `ArgumentParser` object.

    Returns:
        argparse.ArgumentParser: The configured `argparse` parser object for the
        'min-box' subcommand.

    The 'min-box' subcommand accepts the following arguments:

    positional arguments:
        pdb_path    Path to the input PDB file whose box size will be minimized.
                    This argument is optional. If not provided, the program might
                    expect the path through other configuration mechanisms.

    optional arguments:
        --output    Path to the new PDB file where the minimized box structure
                    will be saved. If not provided, the output might overwrite
                    the input file or be handled in another way by the underlying
                    function.

    The function sets the default action for this subparser to be the
    `cli_minimize_box` function, which will be called when the user invokes
    the 'min-box' subcommand.
    """
    parser = subparsers.add_parser(
        "min-box", description="Minimize box size by rotating system."
    )
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
    parser.set_defaults(func=lambda args: cli_minimize_box(args, parser))
    return parser


def cli_minimize_box(args: argparse.Namespace, parser: argparse.ArgumentParser):
    r"""Command-line interface function to rotate a protein to minimize its box volume.

    This function serves as the entry point when the user executes the 'min-box'
    subcommand of the Simlify CLI. It receives the parsed command-line arguments
    and the argument parser object. It extracts the paths to the input PDB file
    and the desired output PDB file and then calls the `run_minimize_box` function
    from the `simlify.structure.pdb.orientation` module to perform the box
    minimization.

    Args:
        args: An `argparse.Namespace` object containing the parsed command-line
            arguments for the 'min-box' subcommand.
        parser (argparse.ArgumentParser): The argument parser object for the
            'min-box' subcommand, used to display help messages if necessary (though
            not explicitly used in this function).

    Returns:
        None

    The function retrieves the following arguments from the `args` object:
        - `pdb_path`: Path to the input PDB file.
        - `output`: Path to the output PDB file.

    It then directly calls the `run_minimize_box` function with these paths to
    perform the box minimization operation.
    """
    run_minimize_box(args.pdb_path, args.output)
