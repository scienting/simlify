from typing import Any

import argparse
import sys

import yaml
from loguru import logger

from simlify import enable_logging


def setup_logging(args: argparse.Namespace) -> None:
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
