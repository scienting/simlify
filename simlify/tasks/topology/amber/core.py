# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

"""
Provides a class for generating Amber topology and coordinate files using the `tleap` program.

This module defines the `AmberTopoGen` class, which inherits from
`simlify.simulation.topo.TopoGen`. It implements the logic to create input files for
`tleap`, execute `tleap`, and parse its output to generate Amber topology (`.prmtop`)
and coordinate (`.inpcrd`) files. It also handles the loading of force fields,
solvation, and ion addition based on the provided `SimlifyConfig`.
"""

from typing import Any

import os
import subprocess
import tempfile
from collections.abc import Iterable

from loguru import logger

from simlify.configs import SimlifyConfig
from simlify.simulation.topo import TopoGen
from simlify.structure.solvent import get_ion_counts
from simlify.utils import simple_generator

FF_WATER_SOLVENT_BOX_MAP: dict[str, Any] = {
    "tip3p": "TIP3PBOX",
    "tip4p": "TIP4PBOX",
    "tip4pew": "TIP4PEWBOX",
    "tip5p": "TIP5PBOX",
    "opc": "OPCBOX",
    "opc3": "OPC3BOX",
    "pol3": "POL3BOX",
    "spce": "SPCBOX",
}
r"""Maps `ff_water` in simulation contexts to tleap box types."""

TLEAP_PATH = os.environ.get("TLEAP_PATH", "tleap")
r"""Path to tleap executable.

You can specify this by setting the path to the `TLEAP_PATH` environmental variable.
For example:

```bash
export TLEAP_PATH="~/miniconda3/envs/metalflare-dev/bin/tleap
```
"""


class AmberTopoGen(TopoGen):
    r"""Standardized framework for generating topology files for Amber
    simulations using tleap."""

    def __init__(self):
        pass

    @classmethod
    def ff_lines(cls, simlify_config: SimlifyConfig) -> Iterable[str]:
        """Prepares `tleap` commands to load specified force fields.

        This method generates a list of `source leaprc.` commands for `tleap`
        based on the force field settings in the provided `SimlifyConfig`.
        It checks for the presence of protein, water, DNA, RNA, GLYCAM, lipid,
        small molecule, and ion force fields in the configuration and creates
        the corresponding `source` or `loadAmberParams` commands.

        Args:
            simlify_config: A `SimlifyConfig` object containing the simulation
                context and force field settings.

        Returns:
             list of strings, where each string is a `tleap`
            command to load a specific force field.

        Examples:
            >>> from simlify.configs import SimlifyConfig
            >>> config = SimlifyConfig(
            ...     engine={
            ...         "ff": {
            ...             "protein": "ff14SB",
            ...             "water": "tip3p",
            ...             "small_molecule": "gaff2",
            ...         }
            ...     }
            ... )
            >>> list(AmberTopoGen.ff_lines(config))
            ['source leaprc.protein.ff14SB', 'source leaprc.water.tip3p', 'source leaprc.gaff2']
        """
        logger.info("Creating tleap commands for loading force fields")
        tleap_lines = []
        if isinstance(simlify_config.engine.ff.protein, str):
            tleap_lines.append(
                f"source leaprc.protein.{simlify_config.engine.ff.protein}"
            )
        if isinstance(simlify_config.engine.ff.water, str):
            tleap_lines.append(f"source leaprc.water.{simlify_config.engine.ff.water}")
        if isinstance(simlify_config.engine.ff.dna, str):
            tleap_lines.append(f"source leaprc.DNA.{simlify_config.engine.ff.dna}")
        if isinstance(simlify_config.engine.ff.rna, str):
            tleap_lines.append(f"source leaprc.RNA.{simlify_config.engine.ff.rna}")
        if isinstance(simlify_config.engine.ff.glycam, str):
            tleap_lines.append(f"source leaprc.{simlify_config.engine.ff.glycam}")
        if isinstance(simlify_config.engine.ff.lipid, str):
            tleap_lines.append(f"source leaprc.{simlify_config.engine.ff.lipid}")
        if isinstance(simlify_config.engine.ff.small_molecule, str):
            tleap_lines.append(
                f"source leaprc.{simlify_config.engine.ff.small_molecule}"
            )
        if isinstance(simlify_config.engine.ff.ions, str):
            tleap_lines.append(
                f"loadAmberParams frcmod.{simlify_config.engine.ff.ions}"
            )
        return tleap_lines

    @classmethod
    def _parse_logs(cls, log_lines: Iterable[str]) -> dict[str, Any]:
        """Parses information from the `tleap.log` file.

        This method iterates through the lines of the `tleap.log` file and extracts
        relevant information such as duplicate atoms, unknown residues, box dimensions,
        volume, mass, density, the number of solvent molecules added, and the net charge
        of the system.

        Args:
            log_lines: An iterable of strings, where each string represents a line
                from the `tleap.log` file.

        Returns:
            A dictionary containing the parsed information from the
                log file. The dictionary may contain the following keys:

                -   `duplicate_atoms`: A list of dictionaries, where each dictionary
                    contains the `residue_number` and `atom_type` of duplicate atoms.
                -   `unknown_residues`: A list of dictionaries, where each dictionary
                    contains the `residue_name`, `residue_number`, and `residue_type`
                    of unknown residues.
                -   `box_volume`: The total volume of the simulation box in Angstroms
                    cubed.
                -   `box_mass`: The total mass of the system in atomic mass units (amu).
                -   `box_density`: The density of the system in grams per cubic
                    centimeter (g/cc).
                -   `solvent_molecules_num`: The number of solvent molecules added
                    to the system.
                -   `charge_net`: The net charge of the system.
        """
        logger.info("Parsing tleap log")
        # Initializing information that is not always in tleap.log
        tleap_info: dict[str, Any] = {"duplicate_atoms": [], "unknown_residues": []}
        tleap_generator = simple_generator(log_lines)
        for line in tleap_generator:
            # -- residue 23: duplicate [ ND2] atoms (total 2)
            if "-- residue " in line and "duplicate" in line:
                residue_number = int(line.split()[2][:-1])
                atom_type = line.split("]")[0].split("[")[-1].strip()
                tleap_info["duplicate_atoms"].append(
                    {"residue_number": residue_number, "atom_type": atom_type}
                )
            # Unknown residue: CRO   number: 64   type: Nonterminal
            if "Unknown residue:" in line:
                line_split = line.split()
                tleap_info["unknown_residues"].append(
                    {
                        "residue_name": line_split[2],
                        "residue_number": line_split[4],
                        "residue_type": line_split[-1],
                    }
                )
            # Total vdw box size:                   66.262 69.436 79.756 angstroms.
            # Volume: 366956.734 A^3
            # Mass > 183192.002 amu,  Density > 0.829 g/cc
            # Added 8761 residues.
            if "Total vdw box size: " in line:
                line = next(tleap_generator)
                tleap_info["box_volume"] = float(line.split()[1])
                line = next(tleap_generator)
                tleap_info["box_mass"] = float(line.split()[2])
                tleap_info["box_density"] = float(line.split()[-2])
                while "Added" not in line:
                    line = next(tleap_generator)
                tleap_info["solvent_molecules_num"] = int(line.split()[1])

            # Total unperturbed charge:  -7.000000
            # Total perturbed charge:    -7.000000
            if "Total unperturbed charge:" in line:
                tleap_info["charge_net"] = float(line.strip().split()[-1])

        logger.debug("Parsed tleap information: {}", tleap_info)
        return tleap_info

    @classmethod
    def _start_tleap_lines(
        cls, path_structure: str, simlify_config: SimlifyConfig
    ) -> list[str]:
        """Generates the initial lines for the `tleap` input file.

        This method prepares the first part of the `tleap` input script, which
        includes loading the specified force fields and any additional lines
        provided in the `simlify_config.topology.append_lines`. It then loads
        the molecular structure from the given `path_structure`.

        Args:
            path_structure: The path to the molecular structure file (e.g., PDB).
            simlify_config: The `SimlifyConfig` object containing simulation settings.

        Returns:
            A list of strings, where each string is a line for the
                `tleap` input file.
        """
        # Prepare tleap_lines
        tleap_lines: list[str] = []
        tleap_lines.extend(cls.ff_lines(simlify_config))
        tleap_lines.extend(simlify_config.topology.append_lines)

        tleap_lines.append(f"mol = loadpdb {path_structure}")
        return tleap_lines

    @staticmethod
    def _write_input(tleap_lines, tmp_input):
        """Writes the generated `tleap` input lines to a temporary file.

        Args:
            tleap_lines: A list of strings, where each string is a line for the
                `tleap` input file.
            tmp_input: A temporary file object opened in write mode.
        """
        logger.debug("Writing tleap input to {}", tmp_input.name)
        tleap_string = "\n".join(tleap_lines)
        tmp_input.write(tleap_string)
        logger.debug("tleap input:\n{}", tleap_string)
        tmp_input.close()  # Need to close before running the command.

    @staticmethod
    def _run_tleap(tmp_input, simlify_config):
        """Executes the `tleap` program using the generated input file.

        Args:
            tmp_input: A temporary file object containing the `tleap` input.
            simlify_config: The `SimlifyConfig` object containing simulation settings,
                used to determine the working directory for `tleap`.

        Returns:
            An object containing information about the completed `tleap` process,
                including its return code and output.
        """
        logger.info("Running tleap")
        tleap_command = [TLEAP_PATH, "-f", tmp_input.name]
        logger.debug("tleap command: {}", tleap_command)
        completed_process = subprocess.run(
            tleap_command,
            capture_output=True,
            text=True,
            check=False,
            cwd=simlify_config.run.dir_work,
        )
        os.remove(tmp_input.name)  # Remove temporary input file.
        logger.debug("tleap output:\n{}", completed_process.stdout)
        if completed_process.returncode != 0:
            logger.error("tleap failed!")
            logger.debug("tleap errors:\n{}", completed_process.stderr)
        return completed_process

    @classmethod
    def dry_run(
        cls,
        path_structure: str,
        simlify_config: SimlifyConfig,
        **kwargs: Any,
    ) -> dict[str, Any]:
        """Performs a dry run of `tleap` to obtain preliminary information about the system.

        This method executes `tleap` with a minimal input script that loads the
        structure, solvates it, saves a temporary PDB file, and calculates the
        total charge. The output log is then parsed to extract information such as
        box dimensions, volume, mass, density, the number of solvent molecules, and
        the net charge. This information is then used to determine the number of ions
        needed to neutralize the system or achieve a desired ionic concentration.

        Args:
            path_structure: Path to the structure file (must be a PDB file for Amber).
            simlify_config: Simlify configuration object.

        Returns:
            A dictionary containing keyword arguments to be passed
                into the `run` method. This dictionary includes information about the
                number of anions and cations required for the system.

        Examples:
            The base `tleap` input file generated for a dry run with standard
            Amber protein context might look like this:

            ```
            source leaprc.protein.ff19SB
            source leaprc.water.opc3
            # User-defined lines from simlify_config.topology.append_lines would be here
            mol = loadpdb <path_structure>
            solvatebox mol OPC3BOX 10.0
            savepdb mol <temporary_file>.pdb
            charge mol
            quit
            ```
        """
        tmp_pdb = tempfile.NamedTemporaryFile(mode="r", suffix=".pdb", delete=True)
        tmp_input = tempfile.NamedTemporaryFile(mode="w+", suffix=".in", delete=False)

        tleap_lines = cls._start_tleap_lines(path_structure, simlify_config)
        solv_box_name = FF_WATER_SOLVENT_BOX_MAP[simlify_config.engine.ff.water]
        tleap_lines.append(
            f"solvatebox mol {solv_box_name} {str(simlify_config.solution.solvent_padding)}"
        )
        tleap_lines.append(f"savepdb mol {tmp_pdb.name}")
        logger.debug("Setting PDB output to {}", tmp_pdb.name)
        tleap_lines.extend(["charge mol", "quit"])

        cls._write_input(tleap_lines, tmp_input)

        completed_process = cls._run_tleap(tmp_input, simlify_config)

        tleap_info = cls._parse_logs(completed_process.stdout.split("\n"))

        if (
            "charge_net"
            not in tleap_info.keys() | "solvent_molecules_num"
            not in tleap_info.keys()
        ):
            logger.error("tleap could not determine relevant system information")
            logger.error(
                "Please check the tleap output for more information: {}", tmp_input.name
            )
        ion_counts = get_ion_counts(
            simlify_config=simlify_config,
            charge_net=tleap_info["charge_net"],
            n_waters=tleap_info["solvent_molecules_num"],
        )

        tleap_info = tleap_info | ion_counts

        return tleap_info

    @classmethod
    def run(
        cls,
        path_structure: str,
        simlify_config: SimlifyConfig,
        charge_anion_num: int = 0,
        charge_cation_num: int = 0,
        path_tleap_pdb: None | str = None,
        **kwargs: Any,
    ) -> dict[str, Any]:
        """Runs `tleap` to generate Amber topology and coordinate files for the system.

        This method takes the path to the structure file, the `SimlifyConfig` object,
        and keyword arguments (typically obtained from the `dry_run` method) to
        generate the final `tleap` input script. It adds ions to neutralize the system
        or achieve a desired concentration, solvates the system, and then saves the
        Amber topology (`.prmtop`) and coordinate (`.inpcrd`) files to the specified
        output paths.

        Args:
            path_structure: Path to the structure file (must be a PDB file for Amber).
            simlify_config: Simlify configuration object.
            charge_anion_num: The number of anions to add.
            charge_cation_num: The number of cations to add.
            path_tleap_pdb: Specify the path to the prepared system as a
                PDB file after `tleap` is finished. If not provided, a
                temporary file is used.

        Returns:
            A dictionary containing information parsed from the
                `tleap` log file.
        """
        # Handle file paths
        if not path_tleap_pdb:
            tmp_tleap_pdb = tempfile.NamedTemporaryFile(
                mode="r", suffix=".pdb", delete=True
            )
            path_tleap_pdb = tmp_tleap_pdb.name
        logger.debug("Setting PDB output to {}", path_tleap_pdb)
        path_topo_write = os.path.join(
            simlify_config.run.dir_work,
            simlify_config.run.dir_input,
            simlify_config.engine.cli.prmtop,
        )
        path_coord_write = os.path.join(
            simlify_config.run.dir_work,
            simlify_config.run.dir_input,
            simlify_config.engine.cli.inpcrd,
        )
        tmp_input = tempfile.NamedTemporaryFile(mode="w+", suffix=".in", delete=False)

        # Prepare and write tleap input
        tleap_lines = cls._start_tleap_lines(path_structure, simlify_config)
        tleap_lines.append(
            f"addIons2 mol {simlify_config.solution.charge_anion_identity} {charge_anion_num}"
        )
        tleap_lines.append(
            f"addIons2 mol {simlify_config.solution.charge_cation_identity} {charge_cation_num}"
        )
        solv_box_name = FF_WATER_SOLVENT_BOX_MAP[simlify_config.engine.ff.water]
        tleap_lines.append(
            f"solvatebox mol {solv_box_name} {str(simlify_config.solution.solvent_padding)}"
        )
        tleap_lines.append(f"savepdb mol {path_tleap_pdb}")
        tleap_lines.append(f"saveamberparm mol {path_topo_write} {path_coord_write}")
        tleap_lines.extend(["charge mol", "quit"])
        cls._write_input(tleap_lines, tmp_input)

        completed_process = cls._run_tleap(tmp_input, simlify_config)
        tleap_info = cls._parse_logs(completed_process.stdout.split("\n"))
        return tleap_info
