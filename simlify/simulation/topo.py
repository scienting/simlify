# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

"""Module defining an abstract base class for topology generators and a function to run topology generation."""

from typing import Any

from abc import ABC, abstractmethod

from simlify import SimlifyConfig
from simlify.utils import get_obj_from_string


class TopoGen(ABC):
    r"""Abstract base class defining the standardized interface for topology file
    generators.

    This class serves as a blueprint for creating topology generators for various molecular
    simulation packages. Subclasses of `TopoGen` are expected to implement the `dry_run`
    and `run` methods, which handle the preliminary information gathering and the actual
    topology file generation, respectively. This ensures a consistent and extensible
    framework for managing different topology generation workflows within Simlify.

    Subclasses should be designed to handle the specific requirements and file formats
    of the target simulation package (e.g., Amber, GROMACS).
    """

    def __init__(self):
        pass

    @classmethod
    @abstractmethod
    def dry_run(
        cls,
        path_structure: str,
        simlify_config: SimlifyConfig,
        **kwargs: Any,
    ) -> dict[str, Any]:
        """Perform a preliminary dry run to gather necessary information before
        topology generation.

        This abstract class method is intended to be implemented by subclasses.
        It should perform a lightweight analysis of the input structure file and the
        Simlify configuration to determine any parameters or information required for
        the full topology generation process. This step can include tasks like
        identifying residue types, counting atoms, or checking for specific force
        field requirements.

        Args:
            path_structure: The file path to the molecular structure file
                (e.g., PDB, GRO) for which the topology will be generated.
            simlify_config: An instance of the Simlify configuration object,
                providing access to global and simulation-specific settings.

        Returns:
            A dictionary containing keyword arguments that should be passed to the
                `run` method of the same topology generator class. This dictionary can
                include any information gathered during the dry run that is needed for
                the actual topology file generation.

        Raises:
            NotImplementedError: If this abstract method is called directly. Subclasses must
                provide their own implementation.
        """
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def run(
        cls,
        path_structure: str,
        simlify_config: SimlifyConfig,
        **kwargs: Any,
    ) -> dict[str, Any]:
        """Generate the topology file based on the provided structure and configuration.

        This abstract class method is intended to be implemented by subclasses. It should
        perform the main task of generating the topology file for the specified simulation
        package. This typically involves reading the structure file, applying force field
        parameters, defining bonds, angles, dihedrals, and other interactions, and writing
        the output topology file in the appropriate format.

        Args:
            path_structure: The file path to the molecular structure file.
            simlify_config: An instance of the Simlify configuration object.

        Returns:
            A dictionary containing information about the generated topology,
                such as the path to the generated file or any other relevant metadata.
                The contents of this dictionary will depend on the specific topology
                generator.

        Raises:
            NotImplementedError: If this abstract method is called directly.
                Subclasses must provide their own implementation.
        """
        raise NotImplementedError


def run_gen_topo(
    path_structure: str,
    import_string: str,
    simlify_config: SimlifyConfig,
) -> dict[str, Any]:
    """Driver function to instantiate and run a topology generator.

    This function takes the path to a structure file, an import string specifying the
    topology generator class to use, and the Simlify configuration. It dynamically
    imports the specified topology generator class, performs a dry run to gather
    preliminary information, and then executes the full topology generation.

    Args:
        path_structure: The file path to the molecular structure file for which
            the topology will be generated.
        import_string: A string representing the import path to the topology
            generator class. This string should be in the format that can be used by
            `simlify.utils.get_obj_from_string` (e.g., `"module.submodule.ClassName"`).
            For example, to use the Amber topology generator, the import string might be
            `"simlify.simulation.amber.topo.AmberTopoGen"`.
        simlify_config: An instance of the Simlify configuration object.

    Returns:
        A dictionary containing information about the generated topology,
            as returned by the `run` method of the instantiated topology generator class.
            This typically includes the path to the generated topology file and other
            relevant metadata.

    Raises:
        ImportError: If the `import_string` does not point to a valid class.
        AttributeError: If the imported class does not have the expected `dry_run` and
            `run` methods.
        Exception: Any exception raised by the `dry_run` or `run` methods of the
            topology generator class.

    Examples:
        To generate an Amber topology for the structure file "input.pdb" using the
        `AmberTopoGen` class:

        >>> from simlify import SimlifyConfig
        >>> config = SimlifyConfig()
        >>> topo_info = run_gen_topo(
        ...     path_structure="input.pdb",
        ...     import_string="simlify.simulation.amber.topo.AmberTopoGen",
        ...     simlify_config=config,
        ... )
        >>> print(topo_info)
        {'topology_file': '/path/to/input.prmtop', 'coordinate_file': '/path/to/input.inpcrd'}
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
