from typing import Any

import argparse
import sys

import yaml
from loguru import logger

from simlify import __version__, enable_logging


def setup_logging(args):
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
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            return dict(yaml.safe_load(f))
    except Exception as e:
        logger.error(f"Error loading YAML file: {e}")
        sys.exit(1)


def cli_main():
    parser = argparse.ArgumentParser(description="Simlify CLI")
    parser.add_argument("-v", action="store_true", help="More log verbosity.")
    parser.add_argument("-vv", action="store_true", help="Even more log verbosity.")
    parser.add_argument("--logfile", help="Specify a file to write logs to.")
    parser.add_argument("--config", help="Path to YAML configuration file.")

    subparsers = parser.add_subparsers(dest="command", help="Available commands:")

    args = parser.parse_args()
    setup_logging(args)

    logger.info(f"Simlify v{__version__} by OASCI <us@oasci.org>")
