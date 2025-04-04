from typing import Any

import argparse
from abc import ABC, abstractmethod

from simlify import SimlifyConfig
from simlify.utils import get_obj_from_string


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
