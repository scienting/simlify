import os
import sys

from loguru import logger

from simlify import __version__
from simlify.cli.config import load_yaml_config, setup_logging
from simlify.cli.parsers import create_parser, print_help_if_no_subcommand


def cli_main():
    parser = create_parser()
    args = parser.parse_args()

    # Print help if no subcommand provided (top-level or nested)
    print_help_if_no_subcommand(parser, args)

    setup_logging(args)
    logger.info(f"Simlify v{__version__} by OASCI <us@oasci.org>")

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
