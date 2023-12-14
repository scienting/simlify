from typing import Any

import os
import subprocess
import tempfile
from collections.abc import Iterable

from loguru import logger

from ...structure.solvent import get_ion_counts
from ...utils import simple_generator
from ..contexts import SimContextManager
from ..topo import TopoGen

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
    r"""Standardized framework for generating topology files."""

    def __init__(self):
        pass

    @classmethod
    def ff_lines(cls, sim_context_manager: SimContextManager) -> Iterable[str]:
        """Prepare to use force fields in the topology file.

        Args:
            sim_context_manager: A simulation context for system preparation.

        Returns:
            `source leaprc.` commands for tleap.

        Examples:
            ```python
            amber_context = {
                "ff_protein": "ff14SB", "ff_water": "tip3p", "ff_small_molecule": "gaff2"
            }
            sim_context_manager = SimulationContextManager(**amber_context)
            tleap_lines = get_source_ff_lines(sim_context_manager)
            ```

            would result in

            ```text
            ["source leaprc.protein.ff14SB", "source leaprc.water.tip3p", "source leaprc.gaff2"]
            ```
        """
        logger.info("Creating tleap commands for loading force fields")
        context = sim_context_manager.get()
        tleap_lines = []
        if isinstance(context["ff_protein"], str):
            tleap_lines.append(f"source leaprc.protein.{context['ff_protein']}")
        if isinstance(context["ff_water"], str):
            tleap_lines.append(f"source leaprc.water.{context['ff_water']}")
        if isinstance(context["ff_dna"], str):
            tleap_lines.append(f"source leaprc.DNA.{context['ff_dna']}")
        if isinstance(context["ff_rna"], str):
            tleap_lines.append(f"source leaprc.RNA.{context['ff_rna']}")
        if isinstance(context["ff_glycam"], str):
            tleap_lines.append(f"source leaprc.{context['ff_glycam']}")
        if isinstance(context["ff_lipid"], str):
            tleap_lines.append(f"source leaprc.{context['ff_lipid']}")
        if isinstance(context["ff_small_molecule"], str):
            tleap_lines.append(f"source leaprc.{context['ff_small_molecule']}")
        if isinstance(context["ff_ions"], str):
            tleap_lines.append(f"loadAmberParams frcmod.{context['ff_ions']}")
        return tleap_lines

    @classmethod
    def _parse_logs(cls, log_lines: Iterable[str]) -> dict[str, Any]:
        """Parse information from any log files.

        Args:
            log_lines: Lines of the `tleap.log` file.

        Returns:
            Parsed information from the log file.

        **Information:**

        -   **`duplicate_atoms`**`: list[dict[str, str]]`

            Each atom in a residue should have a unique type. This list collects information
            of duplicated atom types.

            ```text
            [
                {'residue_number': 23, 'atom_type': 'ND2'},
                {'residue_number': 23, 'atom_type': 'OD1'},
                {'residue_number': 92, 'atom_type': 'NE2'}
            ]
            ```

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

        return tleap_info

    @classmethod
    def _start_tleap_lines(
        cls, path_structure: str, sim_context_manager: SimContextManager
    ) -> list[str]:
        context = sim_context_manager.get()

        # Prepare tleap_lines
        tleap_lines: list[str] = []
        tleap_lines.extend(cls.ff_lines(sim_context_manager))

        extra_lines_topo_gen = context["extra_lines_topo_gen"]
        if extra_lines_topo_gen is None:
            extra_lines_topo_gen = []
        tleap_lines.extend(extra_lines_topo_gen)

        tleap_lines.append(f"mol = loadpdb {path_structure}")
        return tleap_lines

    @staticmethod
    def _write_input(tleap_lines, tmp_input):
        logger.debug("Writing tleap input to {}", tmp_input.name)
        tleap_string = "\n".join(tleap_lines)
        tmp_input.write(tleap_string)
        logger.debug("tleap input:\n{}", tleap_string)
        tmp_input.close()  # Need to close before running the command.

    @staticmethod
    def _run_tleap(tmp_input, dir_work):
        logger.info("Running tleap")
        tleap_command = [TLEAP_PATH, "-f", tmp_input.name]
        logger.debug("tleap command: {}", tleap_command)
        completed_process = subprocess.run(
            tleap_command, capture_output=True, text=True, check=False, cwd=dir_work
        )
        os.remove(tmp_input.name)  # Remove temporary input file.
        logger.debug("tleap output:\n{}", completed_process.stdout)
        if completed_process.returncode != 0:
            logger.error("tleap failed!")
            logger.debug("tleap errors:\n{}", completed_process.stderr)
        return completed_process

    @classmethod
    def dry_run(  # pylint: disable=too-many-arguments
        cls,
        path_structure: str,
        sim_context_manager: SimContextManager,
        dir_work: str | None = None,
    ) -> dict[str, Any]:
        """Perform a dry run to obtain any preliminary information.

        Args:
            path_structure: Path structure file for topology generation. For Amber,
                this must be a PDB file.
            sim_context_manager: Context manager for simulations.
            dir_work: Working directory to generate topology. Useful for
                specifying relative paths.

        Returns:
            Keyword arguments to be passed into `run`.

        Examples:

        The base `tleap` input file is shown below with
        [`AMBER_PROTEIN_STANDARD_CONTEXT`]
        [simulation.amber.contexts.AMBER_PROTEIN_STANDARD_CONTEXT].

        ```bash
        source leaprc.protein.ff19SB
        source leaprc.water.opc3
        <add_lines>
        mol = loadpdb <path_structure>
        solvatebox mol OPC3BOX 10.0
        savepdb mol <temp file>
        charge mol
        quit
        ```
        """
        if dir_work is None:
            dir_work = os.getcwd()
        context = sim_context_manager.get()

        tmp_pdb = tempfile.NamedTemporaryFile(mode="r", suffix=".pdb", delete=True)
        tmp_input = tempfile.NamedTemporaryFile(mode="w+", suffix=".in", delete=False)

        tleap_lines = cls._start_tleap_lines(path_structure, sim_context_manager)
        solv_box_name = FF_WATER_SOLVENT_BOX_MAP[context["ff_water"]]
        tleap_lines.append(
            f"solvatebox mol {solv_box_name} {str(context['solvent_padding'])}"
        )
        tleap_lines.append(f"savepdb mol {tmp_pdb.name}")
        logger.debug("Setting PDB output to {}", tmp_pdb.name)
        tleap_lines.extend(["charge mol", "quit"])

        cls._write_input(tleap_lines, tmp_input)

        completed_process = cls._run_tleap(tmp_input, dir_work)

        tleap_info = cls._parse_logs(completed_process.stdout.split("\n"))

        ion_counts = get_ion_counts(
            sim_context_manager=sim_context_manager,
            charge_net=tleap_info["charge_net"],
            n_waters=tleap_info["solvent_molecules_num"],
        )

        tleap_info = tleap_info | ion_counts

        return tleap_info

    @classmethod
    def run(  # pylint: disable=too-many-arguments
        cls,
        path_structure: str,
        path_topo_write: str,
        path_coord_write: str,
        sim_context_manager: SimContextManager,
        dir_work: str | None = None,
        **kwargs: dict[str, Any],
    ) -> dict[str, Any]:
        r"""Run tleap preparation of a system.

        Args:
            path_structure: Path structure file for topology generation. For Amber,
                this must be a PDB file.
            path_topo_write: Where to write topology file.
            path_coord_write: Where to write coordinate file.
            sim_context_manager: Context manager for simulations.
            dir_work: Working directory to generate topology. Useful for
                specifying relative paths.
            kwargs:
                Keyword arguments used for this specific package.

                Required

                -   `charge_anion_num`
                -   `charge_cation_num`

                Optional

                -   `path_tleap_pdb`: Specify the prepared system as a PDB file after
                tleap is finished.
        """
        if dir_work is None:
            dir_work = os.getcwd()
        context = sim_context_manager.get()

        if "path_tleap_pdb" not in kwargs:
            tmp_tleap_pdb = tempfile.NamedTemporaryFile(
                mode="r", suffix=".pdb", delete=True
            )
            path_tleap_pdb = tmp_tleap_pdb.name
        tmp_input = tempfile.NamedTemporaryFile(mode="w+", suffix=".in", delete=False)

        tleap_lines = cls._start_tleap_lines(path_structure, sim_context_manager)
        n_anions = kwargs["charge_anion_num"]
        n_cations = kwargs["charge_cation_num"]
        tleap_lines.append(
            f"addIons2 mol {context['charge_anion_identity']} {n_anions}"
        )
        tleap_lines.append(
            f"addIons2 mol {context['charge_cation_identity']} {n_cations}"
        )
        solv_box_name = FF_WATER_SOLVENT_BOX_MAP[context["ff_water"]]
        tleap_lines.append(
            f"solvatebox mol {solv_box_name} {str(context['solvent_padding'])}"
        )
        logger.debug("Setting PDB output to {}", path_tleap_pdb)
        tleap_lines.append(f"savepdb mol {path_tleap_pdb}")
        tleap_lines.append(f"saveamberparm mol {path_topo_write} {path_coord_write}")
        tleap_lines.extend(["charge mol", "quit"])

        cls._write_input(tleap_lines, tmp_input)

        completed_process = cls._run_tleap(tmp_input, dir_work)
        # os.remove(tmp_input.name)  # Remove temporary input file.

        tleap_info = cls._parse_logs(completed_process.stdout.split("\n"))

        return tleap_info
