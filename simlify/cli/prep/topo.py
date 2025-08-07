# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

"""
Command-line interface for preparing simulation topology files.
"""

import argparse

from simlify.configs import SimlifyConfig
from simlify.simulation.topo import run_gen_topo


def add_topo_subparser(subparsers):
    r"""Adds the 'topo' subcommand to the Simlify CLI for preparing simulation topology
    files.

    This function configures an `argparse` subparser named 'topo' which allows
    users to specify the path to a molecular structure file, an import string
    pointing to a specific topology generation class, and optional YAML
    configuration files to customize the topology generation process.

    Args:
        subparsers: An `argparse._SubParsersAction` object where the 'topo'
            subparser will be added. This is typically obtained by calling
            `add_subparsers()` on an `ArgumentParser` object.

    Returns:
        The configured `argparse` parser object for the
            'topo' subcommand.

    The 'topo' subcommand accepts the following arguments:

    -   **`structure`:** Path to the molecular structure file (e.g., PDB, GRO)
        for which the topology will be generated.
    -   **`import_string`:** Import string to a topology generation class. This string
        should specify the full path to a class that inherits from
        [`TopoGen`][simulation.topo.TopoGen] and implements the `dry_run` and `run`
        methods. For example: [`AmberTopoGen`][simulation.amber.topo.AmberTopoGen].
    -   **`--yaml`:** One or more paths to YAML configuration files. These files
        will be loaded in the order provided, with later files overriding settings
        from earlier ones. This allows users to configure various aspects of the
        topology generation process.
    -   **`--work`:** Path to the work directory for preparing the topology
        calculations. This option overrides any work directory specified in the YAML
        configuration files. If not provided, the work directory will be determined by
        the configuration.

    The function sets the default action for this subparser to be the
    `cli_topo` function, which will be called when the user invokes the 'topo'
    subcommand.
    """
    parser = subparsers.add_parser("topo", description="Prepare simulation topology.")
    parser.add_argument(
        "structure",
        type=str,
        nargs="?",
        help="Path to a structure file.",
    )
    parser.add_argument(
        "--import-string",
        type=str,
        nargs="?",
        help="Import string to a topology generation class.",
    )
    parser.add_argument(
        "--yaml",
        type=str,
        nargs="+",
        help="Paths to YAML files to use in decreasing precedence.",
    )
    parser.add_argument(
        "--work",
        type=str,
        nargs="?",
        help="Work directory for preparing calculations. This overrides the YAML files.",
        default=None,
    )
    parser.set_defaults(func=lambda args: cli_topo(args, parser))
    return parser


def cli_topo(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    r"""Command-line interface function to prepare simulation topology files.

    This function serves as the entry point when the user executes the 'topo'
    subcommand of the Simlify CLI. It receives the parsed command-line arguments
    and the argument parser object. It initializes a
    [`SimlifyConfig`][config.SimlifyConfig] object and loads configuration settings
    from any YAML files provided via the `--yaml` argument, with later files taking
    precedence. It then checks if a work directory was specified using the `--work`
    option and overrides the work directory in the
    [`SimlifyConfig`][config.SimlifyConfig] if it was provided.
    Finally, it calls [`run_gen_topo`][simulation.topo.run_gen_topo]
    to perform the topology generation.

    Args:
        args: An `argparse.Namespace` object containing the parsed command-line
            arguments for the 'topo' subcommand.
        parser (argparse.ArgumentParser): The argument parser object for the
            'topo' subcommand, used to display help messages if necessary (though
            not explicitly used in this function).

    The function retrieves the following arguments from the `args` object:

    - **`structure`:** Path to the input structure file.
    - **`import_string`:** Import string for the topology generation class.
    - **`yaml`:** A list of paths to YAML configuration files.
    - **`work`:** Optional path to the work directory.

    It initializes a `SimlifyConfig`, loads YAML configurations, overrides the
    work directory if provided, and then calls `run_gen_topo` with the structure
    path, import string, and the configured `SimlifyConfig` object.
    """
    simlify_config = SimlifyConfig()
    if args.yaml is None:
        args.yaml = []
    for yaml_path in reversed(args.yaml):
        simlify_config.from_yaml(yaml_path)
    if args.work is not None:
        simlify_config.run.dir_work = args.work
    run_gen_topo(
        args.structure,
        args.import_string,
        simlify_config,
    )
