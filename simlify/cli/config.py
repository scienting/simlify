# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

"""
Utility functions for the Simlify command-line interface (CLI), including logging setup and configuration loading.
"""

from typing import Any

import argparse
import sys

import yaml
from loguru import logger

from simlify import enable_logging


def setup_logging(args: argparse.Namespace) -> None:
    """Configures the logging level and output for the Simlify CLI based on
    command-line arguments.

    This function takes the parsed command-line arguments and uses them to determine
    the verbosity level of the logs. It supports three levels of verbosity controlled
    by the `-v` and `-vv` flags: INFO (default), DEBUG (`-v`), and TRACE (`-vv`).
    It also allows the user to specify a log file to write the logs to using the
    `--logfile` argument. The logging format is set to display the log level and
    the message.

    Args:
        args: An object containing the parsed command-line arguments, specifically
            looking for `-v`, `-vv`, and `--logfile`.

    The function performs the following steps:

    1. Determines the log level based on the presence of `-v` and `-vv` flags.
       - If `-vv` is present, the log level is set to TRACE (0).
       - If `-v` is present, the log level is set to DEBUG (10).
       - If neither flag is present, the log level defaults to INFO (20).
    2. Calls the `enable_logging` function from the `simlify` package to configure
       the logging system with the determined log level, whether to use standard
       error for output (set to True), the specified log file (if any), and the
       desired log message format.
    """
    if args.vv:
        log_level = 0  # TRACE
    elif args.v:
        log_level = 10  # DEBUG
    else:
        log_level = 20  # INFO

    enable_logging(
        log_level,
        True,
        args.logfile,
        log_format="<level>{level: <8}</level> | {message}",
    )


def load_yaml_config(file_path: str) -> dict[str, Any]:
    """Loads configuration data from a YAML file.

    This function takes the path to a YAML file as input, attempts to open and
    parse the file using the `yaml` library's safe loader, and returns the
    configuration as a Python dictionary. If any error occurs during the file
    reading or parsing process, an error message is logged using `loguru`, and
    the program exits with an error code.

    Args:
        file_path: The path to the YAML configuration file to be loaded.

    Returns:
        A dictionary representing the configuration loaded from the YAML file. The keys
            of the dictionary are strings, and the values can be of any type supported
            by YAML.

    Raises:
        SystemExit: If there is an error loading the YAML file (e.g., file not
            found, invalid YAML syntax), the function will log an error and
            terminate the program with an exit code of 1.
    """
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            return dict(yaml.safe_load(f))
    except Exception as e:
        logger.error(f"Error loading YAML file: {e}")
        sys.exit(1)
