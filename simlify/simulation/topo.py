from typing import Any

import argparse
from abc import ABC, abstractmethod

from ..utils import get_obj_from_string
from .contexts import SimlifyConfig


class TopoGen(ABC):
    r"""Standardized framework for generating topology files."""

    def __init__(self):
        pass

    @classmethod
    @abstractmethod
    def dry_run(  # pylint: disable=too-many-arguments
        cls,
        path_structure: str,
        simlify_config: SimlifyConfig,
    ) -> dict[str, Any]:
        """Perform a dry run to obtain any preliminary information needed before `run`.

        Args:
            path_structure: Path structure file for topology generation.
            simlify_config: Simlify configuration.

        Returns:
            Keyword arguments to be passed into `run`.
        """
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def run(  # pylint: disable=too-many-arguments
        cls,
        path_structure: str,
        simlify_config: SimlifyConfig,
        **kwargs: dict[str, Any],
    ) -> dict[str, Any]:
        """Generate a topology file.

        Args:
            path_structure: Path structure file for topology generation.
            simlify_config: Simlify configuration.
        """
        raise NotImplementedError


def run_gen_topo(
    path_structure: str,
    import_string: str,
    simlify_config: SimlifyConfig,
) -> dict[str, Any]:
    """Diver function for generating a topology file.

    Args:
        path_structure: Path structure file for topology generation.
        import_string: Import string to a topology generation class. For example,
            [`"simlify.simulation.amber.topo.AmberTopoGen"`]
            [simulation.amber.topo.AmberTopoGen].
        simlify_config: Simlify configuration.
    """

    cls_topo = get_obj_from_string(import_string)

    topo_info_dry_run: dict[str, Any] = cls_topo.dry_run(  # type: ignore
        path_structure,
        simlify_config,
    )
    topo_info: dict[str, Any] = cls_topo.run(  # type: ignore
        path_structure,
        simlify_config,
        **topo_info_dry_run,
    )
    return topo_info


def cli_run_gen_topo():
    r"""Command-line interface for generating a topology file"""
    parser = argparse.ArgumentParser(description="Run tleap")
    parser.add_argument(
        "structure",
        type=str,
        nargs="?",
        help="Path to a structure file.",
    )
    parser.add_argument(
        "import_string",
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
    args = parser.parse_args()
    simlify_config = SimlifyConfig()
    if args.yaml is None:
        args.yaml = []
    for yaml_path in reversed(args.yaml):
        simlify_config.from_yaml(yaml_path)
    if args.work is not None:
        simlify_config.rendering.dir_work = args.work
    run_gen_topo(
        args.structure,
        args.import_string,
        simlify_config,
    )
