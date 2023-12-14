from typing import Any

import argparse
from abc import ABC, abstractmethod

from ..utils import get_obj_from_string
from .contexts import SimContextManager


class TopoGen(ABC):
    r"""Standardized framework for generating topology files."""

    def __init__(self):
        pass

    @classmethod
    @abstractmethod
    def dry_run(  # pylint: disable=too-many-arguments
        cls,
        path_structure: str,
        sim_context_manager: SimContextManager,
        dir_work: str | None = None,
    ) -> dict[str, Any]:
        """Perform a dry run to obtain any preliminary information needed before `run`.

        Args:
            path_structure: Path structure file for topology generation.
            sim_context_manager: Context manager for simulations.
            dir_work: Working directory to generate topology. Useful for
                specifying relative paths.

        Returns:
            Keyword arguments to be passed into `run`.
        """
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def run(  # pylint: disable=too-many-arguments
        cls,
        path_structure: str,
        path_topo_write: str,
        path_coord_write: str,
        sim_context_manager: SimContextManager,
        dir_work: str | None = None,
        **kwargs: dict[str, Any],
    ) -> dict[str, Any]:
        """Generate a topology file.

        Args:
            path_structure: Path structure file for topology generation.
            path_topo_write: Where to write topology file.
            path_coord_write: Where to write coordinate file.
            sim_context_manager: Context manager for simulations.
            dir_work: Working directory to generate topology. Useful for
                specifying relative paths.
        """
        raise NotImplementedError


# pylint: disable-next=too-many-arguments
def run_gen_topo(
    path_structure: str,
    path_topo_write: str,
    path_coord_write: str,
    import_string: str,
    sim_context_manager: SimContextManager,
    dir_work: str | None = None,
) -> dict[str, Any]:
    r"""Diver function for generating a topology file.

    Args:
        path_structure: Path structure file for topology generation.
        path_topo_write: Where to write topology file.
        path_coord_write: Where to write coordinate file.
        import_string: Import string to a topology generation class. For example,
            [`"simlify.simulation.amber.topo.AmberTopoGen"`]
            [simulation.amber.topo.AmberTopoGen].
        sim_context_manager: Context manager for simulations.
        dir_work: Working directory to generate topology. Useful for
            specifying relative paths.
    """
    sim_context_manager.dir_work = dir_work
    cls_topo = get_obj_from_string(import_string)

    topo_info_dry_run: dict[str, Any] = cls_topo.dry_run(  # type: ignore
        path_structure,
        sim_context_manager,
        dir_work,
    )
    topo_info: dict[str, Any] = cls_topo.run(  # type: ignore
        path_structure,
        path_topo_write,
        path_coord_write,
        sim_context_manager,
        dir_work,
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
        "topology",
        type=str,
        nargs="?",
        help="Where to save topology file.",
    )
    parser.add_argument(
        "coordinate",
        type=str,
        nargs="?",
        help="Where to save coordinate file.",
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
        help="Work directory for preparing calculations.",
        default=None,
    )
    args = parser.parse_args()
    sim_context_manager = SimContextManager()
    if args.yaml is None:
        args.yaml = []
    for yaml_path in reversed(args.yaml):
        sim_context_manager.from_yaml(yaml_path)
    run_gen_topo(
        args.structure,
        args.topology,
        args.coordinate,
        args.import_string,
        sim_context_manager,
        args.work,
    )
