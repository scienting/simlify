# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

"""
Command-line interface for preparing SLURM submission scripts for simulations.
"""

import argparse

from simlify import SimlifyConfig
from simlify.simulation.slurm import run_slurm_prep


def add_slurm_subparser(subparsers):
    r"""Adds the 'slurm' subcommand to the Simlify CLI for rendering and writing
    SLURM files.

    This function configures an `argparse` subparser named 'slurm' which allows
    users to specify the directory where the SLURM script should be written, the
    path for the SLURM submission script relative to that directory, and optional
    YAML configuration files to customize the SLURM settings.

    Args:
        subparsers: An `argparse._SubParsersAction` object where the 'slurm'
            subparser will be added. This is typically obtained by calling
            `add_subparsers()` on an `ArgumentParser` object.

    Returns:
        The configured `argparse` parser object for the 'slurm' subcommand.

    The 'slurm' subcommand accepts the following arguments:

    -   **`dir_write`:** Directory to write the SLURM submission script and
        potentially other input files.
    -   **`path_slurm_write`:** Local path to write the SLURM submission script,
        relative to the `dir_write`. For example, if `dir_write` is
        `/home/user/simulations` and `path_slurm_write` is `submit.sh`, the script
        will be written to `/home/user/simulations/submit.sh`.
    -   **`--yaml`:** One or more paths to YAML configuration files. These files will
        be loaded in the order provided, with later files overriding settings from
        earlier ones. This allows users to specify default SLURM settings in a base YAML
        file and then override specific parameters for a particular simulation using
        another YAML file.

    The function sets the default action for this subparser to be the
    `cli_slurm_prep` function, which will be called when the user invokes
    the 'slurm' subcommand.
    """
    parser = subparsers.add_parser("slurm", description="Render and write slurm files.")
    parser.add_argument(
        "dir_write",
        type=str,
        nargs="?",
        help="Directory to write input files.",
    )
    parser.add_argument(
        "path_slurm_write",
        type=str,
        nargs="?",
        help="Local path to write a slurm submission script.",
    )
    parser.add_argument(
        "--yaml",
        type=str,
        nargs="*",
        help="Paths to YAML files to use in decreasing precedence.",
    )
    parser.set_defaults(func=lambda args: cli_extract_atoms(args, parser))
    return parser


def cli_extract_atoms(
    args: argparse.Namespace, parser: argparse.ArgumentParser
) -> None:
    r"""Command-line interface function to prepare SLURM submission scripts.

    This function serves as the entry point when the user executes the 'slurm'
    subcommand of the Simlify CLI. It receives the parsed command-line arguments
    and the argument parser object. It initializes a `SimlifyConfig` object and
    loads configuration settings from any YAML files provided via the `--yaml`
    argument, with later files taking precedence. It then sets the SLURM job name
    from the `args.name_job` attribute (which should be added by a higher-level
    parser or configuration). Finally, it calls the `run_slurm_prep` function from
    the `simlify.simulation.slurm` module to create the SLURM submission script.

    Args:
        args: An `argparse.Namespace` object containing the parsed command-line
            arguments for the 'slurm' subcommand.
        parser: The argument parser object for the 'slurm' subcommand, used to
            display help messages if necessary (though not explicitly used in this
            function).

    The function retrieves the following arguments from the `args` object:

    - **`dir_write`:** The directory where the SLURM script will be written.
    - **`path_slurm_write`:** The path to the SLURM script relative to `dir_write`.
    - **`yaml`:** A list of paths to YAML configuration files.
    - **`name_job`:** The name of the SLURM job (expected to be set elsewhere).

    It initializes a `SimlifyConfig`, loads YAML configurations in reverse order
    of precedence, sets the SLURM job name, and then calls `run_slurm_prep` to
    perform the SLURM script preparation.
    """
    simlify_config = SimlifyConfig()
    if args.yaml is None:
        args.yaml = []
    for yaml_path in reversed(args.yaml):
        simlify_config.from_yaml(yaml_path)
    simlify_config.slurm.job_name = args.name_job
    run_slurm_prep(
        args.dir_write,
        args.path_slurm_write,
        simlify_config,
    )
