# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

"""Defines the main command-line interface (CLI) for the Simlify package.

This module provides the entry point for the Simlify command-line tool,
allowing users to interact with Simlify's functionalities through various
subcommands. It handles argument parsing, configuration loading, logging setup,
and dispatches the execution to the appropriate subcommand based on user input.
"""

import os
import sys

from loguru import logger

from simlify import __version__
from simlify.cli.config import load_yaml_config, setup_logging
from simlify.cli.parsers import create_parser, print_help_if_no_subcommand


def cli_main():
    """The main function that powers the Simlify command-line interface.

    This function orchestrates the entire CLI process, from parsing user-provided
    arguments to executing the requested Simlify functionality. It performs the
    following key steps:

    1.  **Argument Parsing**: Utilizes the `create_parser` function to define and parse
        command-line arguments provided by the user. This includes top-level arguments
        and arguments specific to subcommands.

    2.  **Help Display**: Checks if a subcommand was provided. If not, it prints the
        help message for the top-level parser to guide the user on available options
        and subcommands.

    3.  **Logging Setup**: Configures the logging system using the `setup_logging`
        function, which takes the parsed arguments (specifically logging-related
        arguments) to determine the logging level and output format.

    4.  **Version Information**: Logs the Simlify version and the Scienting Studio organization
        information to provide context about the running software.

    5.  **Configuration Loading (Optional)**: Checks if the user provided a path to a
        YAML configuration file using the `--config` argument.
        -   If a configuration file path is provided, it verifies the existence of the
            file. If the file does not exist, a critical error is logged, and the
            program exits.
        -   If the file exists, it loads the configuration using the `load_yaml_config`
            function.
        -   The loaded configuration is then used to update the parsed arguments.
            Command-line arguments take precedence over configuration file settings.
            Only arguments that are not already set via the command line or are set to
            `None` are updated from the configuration file.

    6.  **Subcommand Dispatch**: Checks if a subcommand function (`func`) is associated
        with the parsed arguments. This association is typically established during the
        argument parsing process for each defined subcommand.
        -   If a subcommand function is found, it is called with the parsed arguments as
            input, effectively executing the functionality associated with that subcommand.
        -   If no subcommand function is found (which should ideally be caught by the
            initial help check), it prints the help message for the top-level parser and
            exits with an error code, indicating that the user needs to provide a valid
            subcommand.

    Raises:
        SystemExit: If the provided configuration file does not exist or if no valid
            subcommand is provided by the user.
    """
    parser = create_parser()
    args = parser.parse_args()

    # Print help if no subcommand provided (top-level or nested)
    print_help_if_no_subcommand(parser, args)

    setup_logging(args)
    logger.info(f"Simlify v{__version__} by Scienting Studio <us@scient.ing>")

    # Load configuration file if provided
    if args.config:
        logger.info(f"User provided YAML config file at: `{args.config}`")
        if not os.path.exists(args.config):
            logger.critical("File does not exist. Exiting.")
            sys.exit(1)
        logger.info("Loading YAML file")
        config = load_yaml_config(args.config)
        # Update args with config, prioritizing command-line arguments
        for key, value in config.items():
            if not hasattr(args, key) or getattr(args, key) is None:
                setattr(args, key, value)

    # Dispatch to the subcommand function if available
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    cli_main()
