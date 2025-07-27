# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

from typing import Literal

from atomea.schemas import Render, YamlIO
from loguru import logger
from pydantic import BaseModel

AMBER_CLI_MAPPING = {
    "mdin": "i",
    "mdout": "o",
    "mdinfo": "inf",
    "prmtop": "p",
    "inpcrd": "c",
    "refc": "refc",
    "mtmd": "mtmd",
    "mdcrd": "x",
    "inptraj": "y",
    "mdvel": "v",
    "mdfrc": "frc",
    "mden": "e",
    "restrt": "r",
    "cpin": "cpin",
    "cprestrt": "cprestart",
    "cpout": "cpout",
    "cein": "cein",
    "cerestrt": "cerestrt",
    "ceout": "ceout",
    "evbin": "evbin",
    "suffix": "suffix",
}
"""
Maps keys from `AmberCLIBase` to Amber command-line options.
"""


class AmberCLIBase(BaseModel, YamlIO, Render):
    def render(self, with_newlines: bool = False) -> list[str]:
        """Render the bash command to run an Amber simulation."""
        line: str = self.compute_platform

        for key, value in self.model_dump(exclude_none=True).items():
            if key == "compute_platform":
                continue
            if key[:4] == "dir_":
                # Ignore keys in Render
                continue
            add_option = f" -{AMBER_CLI_MAPPING[key]} {value}"
            logger.debug(f"Adding option: {add_option}")
            line += add_option

        lines = [line]
        if with_newlines:
            lines = [line + "\n" for line in lines]

        return lines

    compute_platform: Literal[
        "sander", "pmemd", "pmemd.MPI", "pmemd.cuda", "pmemd.cuda.MPI"
    ] = "pmemd"

    mdin: str = "md.in"
    """
    **Sander option:** `-i`

    Path to input file for controlling AMBER calculations and operations. Options
    specified in [`AmberInputsBase`]
    [schemas.workflow.amber.inputs.AmberInputsBase] should be in this file.

    Below is a non-working example of what this file should look like.

    ```text

    &cntrl
        imin=0,
        irest=0,
        ntx=1,
        ntmin=1,
        ntave=0,
        ntwv=0,
        ionstepvelocities=0,
        ntwf=0,
        ntwe=0,
        ntwprt=0,
    &end
    ```
    """

    mdout: str = "md.out"
    """
    **Sander option:** `-o`

    Output user readable state info and diagnostics. `-o stdout` will send output to
    stdout (i.e., the terminal) instead of to a file. This stream will contain the
    following information.

    **Characteristics of the Amber release.**

    ```text

            -------------------------------------------------------
            Amber 22 PMEMD                              2022
            -------------------------------------------------------

    | PMEMD implementation of SANDER, Release 22

    |  Compiled date/time: Thu Apr 14 14:06:37 2022
    ```


    Then it will provide any notes and information about the chosen methods and system.

    Afterwords, system information will be printed every
    [`ntpr`][schemas.workflow.amber.inputs.AmberInputsBase.ntpr] steps.

    ```text
    NSTEP =      500   TIME(PS) =    1021.000  TEMP(K) =   300.64  PRESS =     0.0
    Etot   =   -103950.9351  EKtot   =     19726.5940  EPtot      =   -123677.5291
    BOND   =       650.1569  ANGLE   =      1869.0536  DIHED      =      1260.6252
    UB     =         0.0000  IMP     =         0.0000  CMAP       =       155.1250
    1-4 NB =       787.8972  1-4 EEL =     10500.5192  VDWAALS    =     18651.1452
    EELEC  =   -157552.0512  EHBOND  =         0.0000  RESTRAINT  =         0.0000
    EKCMT  =         0.0000  VIRIAL  =         0.0000  VOLUME     =    317418.4248
                                                        Density    =         1.0361
    Ewald error estimate:   0.1408E-04
    ------------------------------------------------------------------------------
    ```

    The meaning of each term is below for convenience.
    You can find more information in the Amber manual.

    -   `NSTEP`: Current step number.
    -   `TIME(PS)`: Cumulative simulation time in picoseconds.
        The previous simulation time is included here if the simulation is a restart.
    -   `TEMP(K)`: Instantaneous system temperature in Kelvin.
    -   `PRESS`: System pressure.
    -   `Etot`: Total energy as the sum of kinetic and potential energies.
    -   `EKtot`: Total kinetic energy.
    -   `EPtot`: Total potential energy.
    -   `BOND`: Sum of bond energies.
    -   `ANGLE`: Sum of angle energies.
    -   `DIHED`: Sum of dihedral energies.
    -   `UB`: Urey-Bradley energy term in the CHARMM force field.
    -   `IMP`: Improper energy term as in the CHARMM force field.
    -   `CMAP`: Correction map for dihedral energies
        [first introduced in CHARMM](https://doi.org/10.1002/jcc.20082) in the
        [extended in Amber's ff19SB](https://doi.org/10.1021/acs.jctc.9b00591)
        (replaces cosine-based $\phi$/$\psi$ dihedral terms in ff14SB).
    -   `1-4 NB`: 1-4 non-bonded (i.e., van der Waal) energy between atoms separated
        by three consecutive bonds. These interactions scale unlike the standard terms,
        so they are treated separately.
    -   `1-4 EEL`: 1-4 electrostatic energy between atoms separated by three consecutive bonds.
        These interactions scale unlike the standard terms, so they are treated separately.
    -   `VDWAALS`: van der Waals energy.
    -   `EELEC`: Electrostatic energy.
    -   `EHBOND`: Depreciated hydrogen bonding contributions.
    -   `RESTRAINT`: Energy of any `ntr` constraints.
    -   `EKCMT`: ?
    -   `VIRIAL`: ?
    -   `VOLUME`: System volume in Å<sup>3</sup>.
    -   `Density`: System density in g/cm<sup>3</sup>.
    """

    mdinfo: str = "md.info"
    """
    **Sander option:** `-inf`

    Path to store the latest
    [`mdout`][schemas.workflow.amber.cli.AmberCLIBase.mdout] and other
    simulation progress information. An example is shown below.

    ```text

    NSTEP =   494000   TIME(PS) =    2008.000  TEMP(K) =   301.01  PRESS =     0.0
    Etot   =   -103688.2116  EKtot   =     19750.7918  EPtot      =   -123439.0034
    BOND   =       691.6081  ANGLE   =      1937.7561  DIHED      =      1279.1292
    UB     =         0.0000  IMP     =         0.0000  CMAP       =       142.8106
    1-4 NB =       795.8815  1-4 EEL =     10392.2564  VDWAALS    =     18105.7369
    EELEC  =   -156784.1822  EHBOND  =         0.0000  RESTRAINT  =         0.0000
    EKCMT  =         0.0000  VIRIAL  =         0.0000  VOLUME     =    317340.9390
                                                        Density    =         1.0364
    Ewald error estimate:   0.8100E-04
    ------------------------------------------------------------------------------
    | Current Timing Info
    | -------------------
    | Total steps:    500000 | Completed:    494000 ( 98.8%) | Remaining:      6000
    |
    | Average timings for last   25500 steps:
    |     Elapsed(s) =      59.89 Per Step(ms) =       2.35
    |         ns/day =      73.58   seconds/ns =    1174.27
    |
    | Average timings for all steps:
    |     Elapsed(s) =    1142.43 Per Step(ms) =       2.31
    |         ns/day =      74.72   seconds/ns =    1156.31
    |
    |
    | Estimated time remaining:      13.9 seconds.
    ------------------------------------------------------------------------------

    ```
    """

    prmtop: str = "mol.prmtop"
    """
    **Sander option:** `-p`

    Path to a file containing input molecular topology, force field,
    periodic box type, atom and residue names.
    """

    inpcrd: str = "mol.inpcrd"
    """
    **Sander option:** `-c`

    Path to file containing input initial coordinates and (optionally) velocities
    and periodic box size.
    """

    refc: str | None = None
    """
    **Sander option:** `-ref`

    input (optional) reference coords for position restraints; also used for targeted MD
    """

    mtmd: str | None = None
    """
    **Sander option:** `-mtmd`

    input (optional) containing list of files and parameters for targeted MD to
    multiple targets
    """

    mdcrd: str | None = None
    """
    **Sander option:** `-x`

    output coordinate sets saved over trajectory
    """

    inptraj: str | None = None
    """
    **Sander option:** `-y`

    input coordinate sets in trajectory format, when imin=5 or 6
    """

    mdvel: str | None = None
    """
    **Sander option:** `-v`

    output velocity sets saved over trajectory
    """

    mdfrc: str | None = None
    """
    **Sander option:** `-frc`

    output force sets saved over trajectory
    """

    mden: str | None = None
    """
    **Sander option:** `-e`

    output extensive energy data over trajectory (not synchronized with mdcrd or mdvel)
    """

    restrt: str | None = None
    """
    **Sander option:** `-r`

    output final coordinates, velocity, and box dimensions if any - for restarting run
    """

    cpin: str | None = None
    """
    **Sander option:** `-cpin`

    input protonation state definitions
    """

    cprestrt: str | None = None
    """
    **Sander option:** `-cprestart`

    protonation state definitions, final protonation states for restart
    (same format as cpin)
    """

    cpout: str | None = None
    """
    **Sander option:** `-cpout`

    output protonation state data saved over trajectory
    """

    cein: str | None = None
    """
    **Sander option:** `-cein`

    input redox state definitions
    """

    cerestrt: str | None = None
    """
    **Sander option:** `-cerestrt`

    redox state definitions, final redox states for restart (same format as cein)
    """

    ceout: str | None = None
    """
    **Sander option:** `-ceout`

    output redox state data saved over trajectory
    """

    evbin: str | None = None
    """
    **Sander option:** `-evbin`

    input input for EVB potentials
    """

    suffix: str | None = None
    """
    **Sander option:** `-suffix`

    output this string will be added to all unspecified output files that are
    printed (for multisander runs, it will append this suffix to all output files)
    """
