# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

from typing import Literal

from loguru import logger
from pydantic import BaseModel, Field

from simlify.configs.utils import Render, YamlIO


class AmberInputsBase(BaseModel, YamlIO, Render):
    def render(self, with_newlines: bool = False) -> list[str]:
        """Render AMBER input file."""
        lines = ["", "&cntrl"]

        for key, value in self.model_dump(exclude_none=True).items():
            if key in (
                "restraintmask",
                "timask1",
                "scmask1",
                "timask2",
                "scmask2",
                "noshakemask",
                "bellymask",
            ):
                if value == "":
                    continue
                value = f'"{value}"'
            line_to_add = f"    {key}={value},"
            logger.debug(f"Adding input line: {line_to_add.strip()}")
            lines.append(line_to_add)

        lines.extend(["&end", ""])

        if with_newlines:
            lines = [line + "\n" for line in lines]

        return lines

    imin: Literal[0, 1, 5, 6, 7] = Field(default=0)
    """Flag for running the energy minimization procedure.

    <hr>

    **`0`**

    Run molecular dynamics without any minimization.

    Perform molecular dynamics (MD) simulation. This mode generates
    configurations by integrating Newtonian equations of motion, allowing the
    system to sample more configurational space and cross small potential energy
    barriers.

    <hr>

    **`1`**

    Perform energy minimization.

    This mode iteratively moves the atoms down the energy gradient to relax the
    structure until a sufficiently low average gradient is obtained. Minimization is
    useful for preparing a system before MD simulations to remove bad contacts and
    high-energy configurations.

    <hr>

    **`5`**

    Read in a trajectory for analysis using the minimization algorithms.

    Although sander will write energy information in the output files
    (using [`ntpr`][configs.workflow.amber.inputs.AmberInputsBase.ntpr]),
    it is often desirable to calculate the energies of a set of structures at a
    later point. In particular, one may wish to post-process a set of structures using
    a different energy function than was used to generate the structures.
    An example of this is MM-PBSA analysis, where the explicit water is removed and
    replaced with a continuum model.

    If [`imin`][configs.workflow.amber.inputs.AmberInputsBase.imin] is set to
    `5`, sander will read a trajectory file (the `inptraj` argument, specified using
    `-y` on the command line), and will perform the functions described in the `mdin`
    file (e.g., an energy minimization) for each of the structures in this file.
    The final structure from each minimization will be written out to the normal
    `mdcrd` file. If you wish to read in a binary (i.e., NetCDF format) trajectory,
    be sure to set [`ioutfm`][configs.workflow.amber.inputs.AmberInputsBase.ioutfm]
    to `1`. Note that this will result in the output trajectory having NetCDF format
    as well.

    For example, when [`imin`][configs.workflow.amber.inputs.AmberInputsBase.imin]
    is `5` and [`maxcyc`][configs.workflow.amber.inputs.AmberInputsBase.maxcyc] is
    `1000`, sander will minimize each structure in the trajectory for 1000 steps and
    write a minimized coordinate set for each frame to the `mdcrd` file.
    If [`maxcyc`][configs.workflow.amber.inputs.AmberInputsBase.maxcyc] is `1`, the
    output file can be used to extract the energies of each of the coordinate
    sets in the `inptraj` file.

    Trajectories containing box coordinates can be post-processed. In order to read
    trajectories with box coordinates,
    [`ntb`][configs.workflow.amber.inputs.AmberInputsBase.ntb] should be greater
    than `0`.

    <hr>

    **`6`**

    Read in a trajectory for analysis using the molecular dynamics driver.

    Like when [`imin`][configs.workflow.amber.inputs.AmberInputsBase.imin] is `5`,
    this option reads a trajectory file for analysis (the `inptraj` argument,
    specified using `-y` on the command line). Instead of minimizing the potential
    energy of each coordinate set, it instead initiates dynamics from each frame as
    if it were read as a restart file without initial velocities. That is, this
    option is equivalent to outputting each frame as a restart file and starting the
    dynamics with [`irest`][configs.workflow.amber.inputs.AmberInputsBase.irest]
    is `0`. If [`nstlim`][configs.workflow.amber.inputs.AmberInputsBase.nstlim] is
    `0`, then this effectively performs a single point energy for each frame.

    <hr>

    **`7`**

    Listen to the selected internet socket and return energies and forces when
    instructed by an external server.

    When this option is set, sander does not perform MD; instead, it listens
    for messages from a server instructing it to compute the potential energy and
    forces of a system. The server IP address and port number are provided as command
    line arguments `-host` and `-port`. The default values are `-host 127.0.0.1` and
    `-port 31415`. The communication pattern follows the protocol implemented in the
    i-PI software. The i-PI program is a molecular dynamics driver used to perform
    classical and centroid path integral molecular dynamics. When i-PI performs
    classical MD, one can instantiate a single sander process to evaluate the
    potential. However, when i-PI is used to perform PIMD, which involves
    calculating potential energies for several “beads” at each time step, multiple
    instances of sander can be launched to simultaneously evaluate the required
    potential energies. The current implementation of the interface is limited to
    simulations in the NVE an NVT ensembles; therefore, one should launch sander with
    a restart file whose unit cell lattice vectors are consistent with the input
    structure supplied to i-PI.
    """

    irest: Literal[0, 1] = Field(default=0)
    """Flag to restart a simulation from a previously saved restart file.
    Restarting a simulation can help in managing large simulations by breaking them
    into smaller, manageable segments and allows for continued simulation from the
    point of interruption without losing progress.

    <hr>

    **`0`**

    Do not restart the simulation; instead, run as a new simulation. Velocities in the
    input coordinate file, if any, will be ignored, and the time step count will be set
    to `0` (unless overridden by
    [`t`][configs.workflow.amber.inputs.AmberInputsBase.t]).

    <hr>

    **`1`**

    Restart the simulation, reading coordinates and velocities from a previously saved
    restart file. The velocity information is necessary when restarting, so
    [`ntx`][configs.workflow.amber.inputs.AmberInputsBase.ntx] must be
    `5` (for Amber versions much older than 20,
    [`ntx`][configs.workflow.amber.inputs.AmberInputsBase.ntx] must be greater than
    or equal to `4`),
    if [`irest`][configs.workflow.amber.inputs.AmberInputsBase.irest] is `1`.
    """

    ntx: Literal[1, 5] = Field(default=1)
    """Option to read the initial coordinates, velocities, and box size from the
    inpcrd file.

    <hr>

    **`1`**

    File is read with no initial velocity information. Suitable for
    starting simulations with new systems where velocities are generated based
    on [`tempi`][configs.workflow.amber.inputs.AmberInputsBase.tempi]
    Option `1` must be used when one is starting from minimized or model-built
    coordinates.

    <hr>

    **`5`**

    File is read unformatted with no initial velocity information. This is
    less common and mainly used for specific needs when dealing with unformatted
    coordinate files.

    If an MD `restrt` file is specified for `inpcrd` then option `5` is generally used
    (unless you explicitly wish to ignore the velocities that are present).
    """

    ntmin: Literal[0, 1, 2, 3, 4] = Field(default=1)
    """
    Flag for selecting minimization type. Determines the algorithm used for energy
    minimization.

    <hr>

    **`0`**

    Full conjugate gradient minimization. The first four cycles are steepest
    descent at the start of the run and after every nonbonded pair list update.
    Conjugate gradient is slower than steepest descent when far from a minimum
    but becomes more efficient when close to the minimum.

    <hr>

    **`1`**

    For [`ncyc`][configs.workflow.amber.inputs.AmberInputsBase.ncyc] cycles, the
    steepest descent method is used, then the conjugate gradient is switched on.
    This option combines the robustness of steepest descent with the efficiency of
    conjugate gradient, making it a recommended choice for many scenarios.

    <hr>

    **`2`**

    Only the steepest descent method is used. This algorithm is popular
    because it is robust and easy to implement. It is generally effective for
    systems far from equilibrium.

    <hr>

    **`3`**

    The XMIN method is used. This method leverages advanced optimization
    algorithms for more efficient minimization, especially useful for large or
    complex systems.

    <hr>

    **`4`**

    The LMOD method is used. This approach uses Low-Mode Conformational
    Search combined with XMIN for energy relaxation and minimization, particularly
    effective for exploring conformational space in flexible molecules.
    """

    maxcyc: int = Field(default=9999, ge=1)
    """
    Maximum number of minimization cycles allowed. This parameter sets the upper limit
    on the number of cycles the minimization algorithm will perform.
    Values typically range from `1000` to `50000`. Lower values may be sufficient for
    small systems or those close to their minimum energy state, while larger values
    can help ensure convergence in more complex or strained systems.

    tip:
        -   If the system is large or highly strained, consider increasing the value to
            ensure sufficient minimization.
        -   For smaller or well-equilibrated systems, a lower value might be adequate,
            saving computational resources.
        -   Monitor the energy convergence; if the energy stabilizes before reaching
            `maxcyc`, the minimization can be considered complete.
    """

    ncyc: int = Field(default=10, ge=1)
    """
    If [`ntmin`][configs.workflow.amber.inputs.AmberInputsBase.ntmin] is `1`,
    then the minimization method will be switched from steepest descent
    to conjugate gradient after
    [`ncyc`][configs.workflow.amber.inputs.AmberInputsBase.ncyc] cycles.
    Values typically range from `5` to `100`.

    Lower values will switch to conjugate gradient sooner, which can be more efficient
    for systems nearing their minimum energy state. Higher values will keep the
    minimization in steepest descent mode longer, which can be beneficial for systems
    far from equilibrium. The default value is set to `10`, which provides a balance
    between the robustness of steepest descent and the efficiency of conjugate gradient
    methods.

    tip:
        -   Adjust the value based on the size and complexity of your system. For larger
            or more complex systems, a higher value might be necessary to ensure
            sufficient initial minimization.
        -   For systems that are relatively close to their equilibrium state, a lower
            value may expedite the minimization process by switching to the more
            efficient conjugate gradient method sooner.
        -   Monitor the minimization process; if convergence is slow, consider
            increasing `ncyc` to allow more initial steepest descent cycles.
    """

    ig: int = Field(default=-1)
    """
    The seed for the pseudo-random number generator. This affects the starting velocities
    for MD simulations if
    [`tempi`][configs.workflow.amber.inputs.AmberInputsBase.tempi] is nonzero.
    If this is `-1`, a random seed will be computed based on the current date and time.
    This should almost always be `-1` unless you are reproducing a run.

    tip:
        -   Set to `-1` to ensure that each simulation run starts with a different seed,
            providing varied initial conditions for better sampling.
        -   Use a specific integer value if you need to reproduce an exact simulation
            run for debugging or verification purposes.
        -   Ensure that [`tempi`][configs.workflow.amber.inputs.AmberInputsBase.tempi]
            is set appropriately when using this option to affect starting velocities.
    """

    dt: float = Field(default=0.001, gt=0.0)
    """
    The time step in picoseconds. This parameter defines the interval of time between
    each step in the simulation.

    -   With SHAKE
        ([`ntc`][configs.workflow.amber.inputs.AmberInputsBase.ntc] is `2`):
        The maximum value is `0.002` ps (2 fs). SHAKE constrains bond lengths
        involving hydrogen atoms, allowing for a larger time step.
    -   Without SHAKE: The maximum value is `0.001` ps (1 fs). Without bond constraints,
        a smaller time step is needed to maintain numerical stability.

    Since longer simulations are usually desired, the maximum value is typically used.
    However, values lower than the maximum can be used if necessary for the phenomena of
    interest, such as high-frequency motions or other specific needs.

    The use of Hydrogen Mass Repartitioning (HMR) together with SHAKE, allows the time
    step to be increased in a stable fashion by about a factor of two (up to `0.004`)
    by slowing down the high frequency hydrogen motion in the system.
    To use HMR, the masses in the topology file need to be altered before starting
    the simulation. ParmEd can do this automatically with the HMassRepartition option.

    tip:
        Consider the nature of your simulation; if capturing high-frequency motions
        is critical, a smaller time step might be required.
    """

    nstlim: int = Field(default=1, ge=1)
    """
    Number of MD steps to perform. Multiply
    [`nstlim`][configs.workflow.amber.inputs.AmberInputsBase.nstlim] and
    [`dt`][configs.workflow.amber.inputs.AmberInputsBase.dt] to get the duration
    of the MD simulation in picoseconds.

    tip:
        -   For short equilibration runs or quick tests, use lower values of
            [`nstlim`][configs.workflow.amber.inputs.AmberInputsBase.nstlim], such as
            `5000` to `50000`.
        -   For production runs aimed at sampling equilibrium states or studying
            slower processes, higher values of
            [`nstlim`][configs.workflow.amber.inputs.AmberInputsBase.nstlim] are
            recommended, typically in the range of `1000000` to `50000000`
            (corresponding to 2 ns to 100 ns for `dt = 0.002`).
        -   Adjust [`nstlim`][configs.workflow.amber.inputs.AmberInputsBase.nstlim]
            based on the specific requirements of your study, considering both
            the computational resources available and the timescale of the phenomena
            you are investigating.
    """

    nscm: int = Field(default=1000, ge=1)
    """
    Flag for the removal of translational and rotational center-of-mass motion every
    [`nscm`][configs.workflow.amber.inputs.AmberInputsBase.nscm] steps.
    For periodic simulations
    ([`ntb`][configs.workflow.amber.inputs.AmberInputsBase.ntb] is `1` or `2`),
    only the translational center-of-mass motion is removed.
    Reasonable values are between `100` and `2000`. Lower values ensure
    more frequent removal of center-of-mass motion, which can be beneficial
    for stability in certain simulations.

    tip:
        -   Use lower values (closer to `100`) for systems that may experience
            significant drift or when higher precision is required.
        -   Higher values (up to `2000`) can be used for more stable systems or
            to reduce computational overhead in very long simulations.
        -   Monitor the simulation to ensure that the chosen value effectively
            manages the center-of-mass motion without introducing artifacts or
            instability.
    """

    ntr: Literal[0, 1] = Field(default=0)
    """
    Flag for restraining positions of specified atoms using a harmonic potential.
    Ensure that
    [`restraintmask`][configs.workflow.amber.inputs.AmberInputsBase.restraintmask]
    is properly defined to specify the atoms that require constraints.

    <hr>

    **`0`**

    No constraints. The positions of all atoms are free to move according to the
    simulation dynamics. Use `0` for fully flexible simulations where no positional
    restraints are needed.

    <hr>

    **`1`**

    Constrain atoms specified in
    [`restraintmask`][configs.workflow.amber.inputs.AmberInputsBase.restraintmask].
    This applies a harmonic potential to the atoms defined in
    [`restraintmask`][configs.workflow.amber.inputs.AmberInputsBase.restraintmask],
    effectively fixing their positions relative to the rest of the system.
    Use `1` when specific atoms need to be restrained, such as in cases where you want
    to focus on a particular region of the system while keeping another region fixed
    or minimally perturbed.
    """

    restraint_wt: float = Field(default=4.0, ge=0.0)
    """
    The weight (in kcal mol<sup>-1</sup> Å<sup>-2</sup>) when
    [`ntr`][configs.workflow.amber.inputs.AmberInputsBase.ntr] is `1`. The form of
    the restraint is
    $k (\Delta x)^2$ where $\Delta x$ is the deviation of the atom's coordinate
    from the reference position.

    -   Reasonable values are between `0.1` and `10.0`.
        -   Lower values (`0.1` to `1.0`) allow for more movement, providing a looser
            restraint.
        -   Higher values (`5.0` to `10.0`) significantly restrict movement, enforcing
            a tighter restraint.
    -   The default value is set to `4.0`, which provides a moderate restraint,
        balancing stability and flexibility.
    """

    restraintmask: str = Field(default="")
    """
    Strings that specify the restricted atoms when
    [`ntr`][configs.workflow.amber.inputs.AmberInputsBase.ntr] is `1`. To see what
    atoms will be restrained, you can use
    `ambmask -p mol.prmtop -c mol.inpcrd -out pdb -find "RESTRAINT_STRING"`
    in `ambertools`. Here are some examples of `restraintmask`s and their descriptions.

    -   `"!(@H=)"`: Restrain all atoms except hydrogens.
    -   `"!(:WAT) & !(@H=)"`: Restrain all atoms except hydrogens and water molecules.
    -   `"!(:WAT) & !(@H=) & !(:Na+,Cl-)"`: Same as above, but does not
        restrain Na<sup>+</sup> and Cl<sup>-</sup> ions.
    -   `"!(:WAT) & (@C,CA,N,O,O5',P,O3',C3',C4',C5')"`: Restrains protein, DNA,
        and RNA backbone atoms.

    We must include the `!(:WAT)` to avoid restraining oxygen atoms in the water molecules.
    """

    ntb: Literal[0, 1, 2] = Field(default=1)
    """
    Flag for periodic boundary conditions when computing non-bonded interactions.
    Bonds that cross the boundary are not supported.

    <hr>

    **`0`**

    No periodic boundary conditions. This is suitable for simulations where boundary
    effects are not a concern, such as in isolated systems or gas-phase simulations.

    <hr>

    **`1`**

    Constant volume. This maintains a fixed simulation box size, appropriate for
    systems where volume changes are not expected or desired.

    <hr>

    **`2`**

    Constant pressure. This allows the simulation box to fluctuate in size to
    maintain constant pressure, suitable for more realistic simulations of
    condensed-phase systems where pressure control is needed.
    """

    ntf: Literal[1, 2, 3, 4, 5, 6, 7, 8] = Field(default=1)
    """
    Force evaluation type. This parameter determines which interactions are
    considered during the force calculations.

    <hr>

    **`1`**

    All contributions. This option includes all types of interactions in the
    force evaluation and is required for minimization.

    <hr>

    **`2`**

    Ignore bond interactions involving hydrogens. This option is typically used when
    [`ntc`][configs.workflow.amber.inputs.AmberInputsBase.ntc] is `2`, meaning
    constraints are applied to bonds involving hydrogens (e.g., SHAKE algorithm).

    <hr>

    **`3`**

    All bond interactions are omitted. This option is used when
    [`ntc`][configs.workflow.amber.inputs.AmberInputsBase.ntc] is `3`, which
    implies constraints are applied to all bonds.

    <hr>

    **`4`**

    angle involving H-atoms and all bonds are omitted

    <hr>

    **`5`**

    all bond and angle interactions are omitted

    <hr>

    **`6`**

    dihedrals involving H-atoms and all bonds and all angle interactions are omitted

    <hr>

    **`7`**

    all bond, angle and dihedral interactions are omitted

    <hr>

    **`8`**

    all bond, angle, dihedral and non-bonded interactions are omitted
    """

    ntc: Literal[1, 2, 3] = Field(default=1)
    """
    Flag for SHAKE to perform bond length constraints. The SHAKE option should be used
    for most MD calculations. The size of the MD timestep is determined by the
    fastest motions in the system. SHAKE removes the bond stretching freedom, which
    is the fastest motion, and consequently allows a larger timestep to be used.
    For water models, a special "three-point" algorithm is used. Consequently, to
    employ TIP3P set [`ntf`][configs.workflow.amber.inputs.AmberInputsBase.ntf] and
    [`ntc`][configs.workflow.amber.inputs.AmberInputsBase.ntc] to 2.

    Since SHAKE is an algorithm based on dynamics, the minimizer is not aware of
    what SHAKE is doing; for this reason, minimizations generally should be carried
    out without SHAKE. One exception is short minimizations whose purpose is to remove
    bad contacts before dynamics can begin.

    For parallel versions of sander only intramolecular atoms can be constrained.
    Thus, such atoms must be in the same chain of the originating PDB file.

    <hr>

    **`1`**

    No SHAKE. This option does not apply any bond length constraints and is
    required for minimization.

    <hr>

    **`2`**

    Bonds involving hydrogen are constrained. This is the recommended setting for
    molecular dynamics (MD) simulations as it allows for larger time steps while
    maintaining stability.

    <hr>

    **`3`**

    All bonds are constrained. This setting applies constraints to all bonds,
    providing the most rigid structure.

    !!! warning
        Not available for parallel or qmmm runs in sander.
    """

    cut: float = Field(default=8.0, gt=0.0)
    """
    Specifies the nonbonded cutoff in Angstroms. This parameter sets the distance
    beyond which nonbonded interactions (such as van der Waals and electrostatic
    interactions) are not explicitly calculated.

    For Particle Mesh Ewald (PME), this value limits the direct space sum.
    PME is a method used to compute long-range electrostatic interactions
    efficiently by splitting them into short-range (direct space) and long-range
    (reciprocal space) components. The cutoff distance specifies the range within
    which the direct space interactions are computed explicitly, while interactions
    beyond this distance are handled in reciprocal space.

    -   `8.0` Å is commonly used in many simulations as a balance between accuracy and
        computational efficiency. It ensures that the majority of significant
        interactions are captured while keeping the computational cost manageable.
    -   `10.0` Å can provide improved accuracy, especially in systems where
        long-range interactions  are important. However, it increases the computational
        load.

    When using an implicit solvent model, the nonbonded cutoff should be set to
    `9999.0` Å, effectively making it infinite. This is because implicit solvent models
    treat solvation effects as a continuous medium, and all interactions need to be
    considered without a cutoff.
    """

    ntt: Literal[0, 1, 2, 3, 9, 10, 11] = Field(default=3)
    """Switch for temperature scaling.

    <hr>

    **`0`**

    Constant total energy (NVE).

    <hr>

    **`1`**

    Constant temperature, using the weak-coupling algorithm. A single scaling factor
    is used for all atoms. Note that this algorithm just ensures that the total
    kinetic energy is appropriate for the desired temperature; it does nothing to
    ensure that the temperature is even over all parts of the molecule. Atomic
    collisions will tend to ensure an even temperature distribution, but this is not
    guaranteed, and there are many subtle problems that can arise with weak
    temperature coupling.

    Using [`ntt`][configs.workflow.amber.inputs.AmberInputsBase.ntt] of `1` is
    especially dangerous for generalized Born simulations, where there are no
    collisions with solvent to aid in thermalization.
    Other temperature coupling options (especially
    [`ntt`][configs.workflow.amber.inputs.AmberInputsBase.ntt] of `3`) should be
    used instead.

    <hr>

    **`2`**

    Andersen-like temperature coupling scheme, in which imaginary "collisions"
    are performed with heat bath of temperature
    [`temp0`][configs.workflow.amber.inputs.AmberInputsBase.temp0] every
    [`vrand`][configs.workflow.amber.inputs.AmberInputsBase.vrand] steps.

    note:
        In between these "massive collisions", the dynamics is Newtonian. Hence, time correlation
        functions (etc.) can be computed in these sections, and the results averaged over an initial
        canonical distribution. Note also that too high a collision rate (too small a value of vrand)
        will slow down the speed at which the molecules explore configuration space, whereas too
        low a rate means that the canonical distribution of energies will be sampled slowly.

    <hr>

    **`3`**

    Use Langevin dynamics with the collision frequency
    [`gamma_ln`][configs.workflow.amber.inputs.AmberInputsBase.gamma_ln].
    Since Langevin simulations are highly susceptible to "synchronization"
    artifacts, you should explicitly set
    [`ig`][configs.workflow.amber.inputs.AmberInputsBase.ig] to a different value
    every restart (e.g., `-1`).

    <hr>

    **`9`**

    Optimized Isokinetic Nose-Hoover chain ensemble (OIN).
    Constant temperature simulation utilizing Nose-Hoover chains and an isokinetic
    constraint on the particle and thermostat velocities, implemented for use in
    multiple time-stepping methods, namely for 3D-RISM and RESPA.

    <hr>

    **`10`**

    Stochastic Isokinetic Nose-Hoover RESPA integrator.
    Mainly used for RESPA simulations.

    <hr>

    **`11`**

    Stochastic version of Berendsen thermostat, also known as the
    Bussi thermostat. This thermostat samples canonical distribution by scaling
    all velocities to a random temperature probed from canonical distribution.
    """

    tempi: float = Field(default=100.0, gt=0.0)
    """
    Initialization temperature in Kelvin. This parameter sets the initial temperature
    for the system at the start of the simulation.

    -   If [`ntx`][configs.workflow.amber.inputs.AmberInputsBase.ntx] is `1`, the
        initial velocities of the atoms are assigned from a
        Maxwellian distribution corresponding to the temperature
        [`tempi`][configs.workflow.amber.inputs.AmberInputsBase.tempi]. This is
        typically used to start a new simulation where the initial conditions need to
        be defined.
    -   If [`ntx`][configs.workflow.amber.inputs.AmberInputsBase.ntx] is `5`, this
        parameter has no effect because the velocities are read from the input
        coordinates file, meaning the system continues from a previously
        equilibrated state.

    tip:
        Use a lower initial temperature (e.g., `100.0` K) to gradually equilibrate
        the system and avoid potential instabilities or large initial forces.
    """

    temp0: float = Field(default=300.0, gt=0.0)
    """
    Reference temperature at which the system will be kept in Kelvin. This parameter
    defines the target temperature for the simulation, around which the thermostat will
    regulate the system. It is crucial for simulations where temperature control is
    needed to mimic real-world conditions, such as simulations of biological systems
    or materials at specific temperatures.

    The default value is set to `300.0` K, which corresponds to approximately room
    temperature. This is a common target temperature for many biological simulations.
    Other typical values might include `310.0` K for physiological conditions
    (human body temperature) or other specific temperatures relevant to your study.
    """

    gamma_ln: float = Field(default=2.0, gt=0.0)
    """
    The collision frequency, $\gamma$, when
    [`ntt`][configs.workflow.amber.inputs.AmberInputsBase.ntt] is `3`. This parameter
    is used in Langevin dynamics to control the rate of coupling between the system
    and a heat bath, thereby regulating the temperature.

    The default value is set to `2.0` ps<sup>-1</sup>, which is commonly used in many
    simulations to provide effective temperature control without overly damping the
    system's dynamics. Values typically range from `2.0` to `5.0` ps<sup>-1</sup>.
    While the physical collision frequency is about `50` ps<sup>-1</sup>, using such
    high values can overly dampen the system. Smaller values are preferred to maintain
    more realistic dynamics.
    """

    ntp: Literal[0, 1, 2, 3, 4] = Field(default=0)
    """
    Flag for constant pressure dynamics. This parameter controls how pressure is
    managed during the simulation.

    <hr>

    **`0`**

    No pressure scaling. The system is run at constant volume, and no
    adjustments are made to maintain a specific pressure.

    <hr>

    **`1`**

    Isotropic position scaling. This is the recommended setting for most
    simulations as it scales the simulation box uniformly in all directions to
    maintain constant pressure.

    <hr>

    **`2`**

    Anisotropic pressure scaling. This option allows different scaling factors
    for each dimension and can only be used for orthogonal boxes. It is typically
    used in membrane simulations where different surface tensions exist in the
    $x$, $y$, and $z$ directions. Solutes dissolved in water should not use this
    setting as it can introduce artifacts.

    <hr>

    **`3`**

    Molecular dynamics with semiisotropic pressure scaling: this is only available with
    constant surface tension (`csurften` > `0`) and orthogonal boxes. This links the
    pressure coupling in the two directions tangential to the interface.

    <hr>

    **`4`**

    Molecular dynamics towards a targeted volume. This is not for production but
    for modifying the volume of the system, particularly useful for preparing
    replica-exchange molecular dynamics simulations where the shape of each replica
    needs to be the same. When
    [`ntp`][configs.workflow.amber.inputs.AmberInputsBase.ntp] is `4`, the following
    variables in the “ewald” namelist should be set:

    -   `target_n`: Number of target volume iterations to reach the target volume
        (default 100).
    -   `target_a`, `target_b`, `target_c`: the cell dimension of the target volume.
    """

    barostat: Literal[1, 2] = Field(default=1)
    """
    Flag used to control which barostat is used to regulate the pressure during the
    simulation. The barostat setting determines how the simulation box's volume is
    adjusted to maintain the desired pressure. Choosing the right barostat is important
    for the accuracy and stability of the simulation.

    <hr>

    **`1`**

    Berendsen barostat. This method scales the box dimensions and atomic
    coordinates to achieve the desired pressure. It is simpler but can lead to
    less accurate pressure control and may not generate a true NPT ensemble.

    <hr>

    **`2`**

    Monte Carlo barostat. This method uses Monte Carlo moves to adjust the
    box volume and is generally more accurate for maintaining constant pressure.
    It provides better sampling of the NPT ensemble and is recommended for most
    simulations.
    """

    pres0: float = Field(default=1.0, gt=0.0)
    """
    Reference pressure, in bar, at which the system is maintained. This is the target
    pressure for the barostat. This value is almost always used in simulations to
    mimic standard atmospheric conditions.
    """

    mcbarint: int = Field(default=100, ge=1)
    """
    Number of steps between volume change attempts performed as part of the Monte Carlo
    barostat. This determines how frequently the volume of the simulation box is
    adjusted during pressure regulation.
    """

    comp: float = Field(default=44.6, gt=0.0)
    """
    Compressibility of the system when
    [`ntp`][configs.workflow.amber.inputs.AmberInputsBase.ntp] > `0` in units of
    10<sup>-6</sup> bar<sup>-1</sup>. 44.6` x 10<sup>-6</sup> bar<sup>-1</sup>,
    appropriate for water. This value is used in pressure regulation to account for the
    compressibility of the solvent or system being simulated.
    """

    taup: float = Field(default=1.0, gt=0.0)
    """
    Pressure relaxation time in picoseconds when
    [`ntp`][configs.workflow.amber.inputs.AmberInputsBase.ntp] > `0`.
    Recommended values are between `1.0` and `5.0` ps. This parameter controls how
    quickly the pressure adjusts to the target value. Start with `1.0` ps.
    If your simulations are unstable, consider increasing this value.
    """

    ntxo: Literal[1, 2] = Field(default=2)
    """
    Format of the final coordinates, velocities, and box size (if a constant volume or
    pressure run) written to file `restrt`.

    <hr>

    **`1`**

    ASCII.

    <hr>

    **`2`**

    Binary NetCDF file.
    """

    ntwr: int = Field(default=1000, ge=1)
    """
    Every [`ntwr`][configs.workflow.amber.inputs.AmberInputsBase.ntwr] steps
    during dynamics, the `restrt` file will be written, ensuring that recovery from a
    crash will not be so painful. No matter what the value of
    [`ntwr`][configs.workflow.amber.inputs.AmberInputsBase.ntwr], a `restrt` file
    will be written at the end of the run, i.e., after
    [`nstlim`][configs.workflow.amber.inputs.AmberInputsBase.nstlim] steps
    (for dynamics) or
    [`maxcyc`][configs.workflow.amber.inputs.AmberInputsBase.maxcyc] steps
    (for minimization). If
    [`ntwr`][configs.workflow.amber.inputs.AmberInputsBase.ntwr] < `0`, a unique
    copy of the file, `restrt_<nstep>`, is written every `abs(ntwr)` steps.
    This option is useful if for example one wants to run free energy perturbations
    from multiple starting points or save a series of `restrt` files for minimization.
    """

    ntpr: int = Field(default=1000, ge=1)
    """
    Print energy information every
    [`ntpr`][configs.workflow.amber.inputs.AmberInputsBase.ntpr] steps in a
    human-readable form to files `mdout` and `mdinfo`. `mdinfo` is closed and
    reopened each time, so it always contains the most recent energy and temperature.
    """

    ntwx: int = Field(default=0, ge=1)
    """
    Coordinates are written every
    [`ntwx`][configs.workflow.amber.inputs.AmberInputsBase.ntwx] steps to the
    `mdcrd` file. This parameter controls how often the coordinates are saved,
    which can be used for trajectory analysis. If
    [`ntwx`][configs.workflow.amber.inputs.AmberInputsBase.ntwx] is `0`, no coordinate
    trajectory file will be written.
    """

    ioutfm: Literal[0, 1] = Field(default=1)
    """
    Format of coordinate and velocity trajectory files.

    <hr>

    **`0`**

    ASCII.

    <hr>

    **`1`**

    Binary NetCDF files.
    """

    iwrap: Literal[0, 1] = Field(default=0)
    """
    Flag for wrapping coordinates around the periodic boundary.

    <hr>

    **`0`**

    No wrapping will be performed, in which case it is typical to use cpptraj as a
    post-processing program to translate molecules back to the primary box.

    <hr>

    **`1`**

    The coordinates written to the restart and trajectory files will be "wrapped" into
    a primary box. This means that for each molecule, its periodic image closest
    to the middle of the "primary box" (with $x$ coordinates between $0$ and $a$, $y$
    coordinates between $0$ and $b$, and $z$ coordinates between $0$ and $c$) will
    be the one written to the output file. This often makes the resulting structures
    look better visually, but has no effect on the energy or forces. Performing such
    wrapping, however, can mess up diffusion and other calculations.

    For very long runs, setting
    [`iwrap`][configs.workflow.amber.inputs.AmberInputsBase.iwrap] is `1` may be
    required to keep the coordinate output from overflowing the trajectory and
    restart file formats, especially if trajectories are written in ASCII format
    instead of NetCDF (see also the
    [`ioutfm`][configs.workflow.amber.inputs.AmberInputsBase.ioutfm] option).
    """

    nmropt: Literal[0, 1, 2] = Field(default=0)
    """

    <hr>

    **`0`**

    No nmr-type analysis will be done.

    <hr>

    **`1`**

    NMR restraints and weight changes will be read.

    <hr>

    **`2`**

    NMR restraints, weight changes, NOESY volumes, chemical shifts and residual dipolar
    restraints will be read.
    """

    ntave: int = Field(default=0, ge=0.0)
    """Every [`ntave`][configs.workflow.amber.inputs.AmberInputsBase.ntave] steps of
    dynamics, running averages of average energies and fluctuations over the last
    [`ntave`][configs.workflow.amber.inputs.AmberInputsBase.ntave] steps will be
    printed out. A value of `0` disables this printout. Setting
    [`ntave`][configs.workflow.amber.inputs.AmberInputsBase.ntave] to a value
    1/2 or 1/4 of nstlim provides a simple way to look at convergence during the
    simulation.
    """

    ntwv: int = Field(default=0, ge=0.0)
    """
    Every [`ntwv`][configs.workflow.amber.inputs.AmberInputsBase.ntwv] steps,
    the velocities will be written to the `mdvel` file. If
    [`ntwv`][configs.workflow.amber.inputs.AmberInputsBase.ntwv] is `0`, no
    velocity trajectory file will be written. If
    [`ntwv`][configs.workflow.amber.inputs.AmberInputsBase.ntwv] is `-1`, velocities
    will be written to mdcrd, which then becomes a combined coordinate/velocity
    trajectory file, at the interval defined by ntwx. This option is available
    only for binary NetCDF output
    ([`ioutfm`][configs.workflow.amber.inputs.AmberInputsBase.ioutfm] is `1`).
    Most users will have no need for a velocity trajectory file and so can safely
    leave [`ntwv`][configs.workflow.amber.inputs.AmberInputsBase.ntwv] at the default.
    Note that dumping velocities frequently, like forces or coordinates, will introduce
    potentially significant I/O and communication overhead, hurting both performance
    and parallel scaling.
    """

    ionstepvelocities: Literal[0, 1] = Field(default=0)
    """Controls whether to print the half-step-ahead velocities (**`0`**) or
    on-step velocities (**`1`**). The half-step-ahead velocities can potentially
    be used to restart calculations, but the on-step velocities correspond to
    calculated kinetic energy/temperature.
    """

    ntwf: int = Field(default=0, ge=0)
    """
    Every [`ntwf`][configs.workflow.amber.inputs.AmberInputsBase.ntwf] steps, the
    forces will be written to the mdfrc file. If
    [`ntwf`][configs.workflow.amber.inputs.AmberInputsBase.ntwf] is `0`, no force
    trajectory file will be written. If
    [`ntwf`][configs.workflow.amber.inputs.AmberInputsBase.ntwf] is `-1`, forces
    will be written to the mdcrd, which then becomes a combind coordinate/force
    trajectory file, at the interval defined by
    [`ntwx`][configs.workflow.amber.inputs.AmberInputsBase.ntwx].
    This option is available only for binary NetCDF output
    ([`ioutfm`][configs.workflow.amber.inputs.AmberInputsBase.ioutfm] is `1`).
    Most users will have no need for a force trajectory file and so can safely
    leave [`ntwf`][configs.workflow.amber.inputs.AmberInputsBase.ntwf] at the default.
    Note that dumping forces frequently, like velocities or coordinates, will introduce
    potentially significant I/O and communication overhead, hurting both performance
    and parallel scaling.
    """

    ntwe: int = Field(default=0, ge=0)
    """
    Every [`ntwe`][configs.workflow.amber.inputs.AmberInputsBase.ntwe] steps,
    the energies and temperatures will be written to file `mden` in a compact form. If
    [`ntwe`][configs.workflow.amber.inputs.AmberInputsBase.ntwe] is `0` then no
    `mden` file will be written. Note that energies in the `mden` file are not
    synchronized with coordinates or velocities in the `mdcrd` or `mdvel` file(s).
    Assuming identical [`ntwe`][configs.workflow.amber.inputs.AmberInputsBase.ntwe]
    and [`ntwx`][configs.workflow.amber.inputs.AmberInputsBase.ntwx] values the
    energies are one time step before the coordinates (as well as the velocities
    which are synchronized with the coordinates). Consequently, an `mden` file is
    rarely written.
    """

    ntwprt: int = Field(default=0, ge=0)
    """
    The number of atoms to include in trajectory files (`mdcrd` and `mdvel`).
    This flag can be used to decrease the size of the these files, by including only
    the first part of the system, which is usually of greater interest (for instance,
    one might include only the solute and not the solvent).

    <hr>

    **`0`**

    Include all atoms of the system when writing trajectories.

    **`> 0`**

    <hr>

    Include only atoms 1 to
    [`ntwprt`][configs.workflow.amber.inputs.AmberInputsBase.ntwprt] when writing
    trajectories.
    """

    idecomp: Literal[0, 1, 2, 3, 4] = Field(default=0)
    r"""
    Perform energy decomposition according to a chosen scheme. In former
    distributions this option was only really useful in conjunction with `mm_pbsa`,
    where it is turned on automatically if required. Now, a
    decomposition of $\left\langle \partial V / \partial \gamma \right\rangle$ on a
    per-residue basis in thermodynamic integration (TI) simulations is also possible.

    If energy decomposition is requested, residues may be chosen by the
    RRES and/or LRES card. The RES card is used to select the residues about which
    information is written out.

    <hr>

    **`0`**

    Do not decompose energies.

    <hr>

    **`1`**

    Decompose energies on a per-residue basis; 1-4 EEL + 1-4 VDW are added to internal
    (bond, angle, dihedral) energies.

    <hr>

    **`2`**

    Decompose energies on a per-residue basis; 1-4 EEL + 1-4 VDW are added to EEL and
    VDW.

    <hr>

    **`3`**

    Decompose energies on a pairwise per-residue basis; otherwise equivalent to
    [`idecomp`][configs.workflow.amber.inputs.AmberInputsBase.idecomp] is `1`.
    Not available in TI simulations.

    <hr>

    **`4`**

    Decompose energies on a pairwise per-residue basis; otherwise equivalent to
    [`idecomp`][configs.workflow.amber.inputs.AmberInputsBase.idecomp] is `2`.
    Not available in TI simulations.
    """

    ibelly: Literal[0, 1] = Field(default=0)
    """
    Flag for belly type dynamics. If set to `1`, a subset of the atoms in the system
    will be allowed to move, and the coordinates of the rest will be frozen.
    The moving atoms are specified with
    [`bellymask`][configs.workflow.amber.inputs.AmberInputsBase.bellymask].
    This option is not available when
    `igb` > `0`. When belly
    type dynamics is in use, bonded energy terms, vdW interactions, and direct
    space electrostatic interactions are not calculated for pairs of frozen atoms.
    Note that this does not provide any significant speed advantage. Freezing
    atoms can be useful for some applications but is maintained primarily for
    backwards compatibility with older versions of Amber. Most applications should
    use the [`ntr`][configs.workflow.amber.inputs.AmberInputsBase.ntr] variable
    instead to restrain parts of the system to stay close to some initial configuration.
    """

    bellymask: str = Field(default="")
    """String that specifies the moving atoms when
    [`ibelly`][configs.workflow.amber.inputs.AmberInputsBase.ibelly] is `1`.
    """

    dx0: float = Field(default=0.01, ge=0.0)
    """
    The initial step length. If the initial step length is too big then will give a
    huge energy; however the minimizer is smart enough to adjust itself.
    """

    drms: float = Field(default=0.0001, gt=0)
    """
    The convergence criterion for the energy Derivative: minimization will halt when
    the Root-Mean-Square of the Cartesian elements of the gradient of the energy is
    less than this.
    """

    t: float = Field(default=0.0, ge=0.0)
    """The time at the start (psec) this is for your own reference and is not critical.
    Start time is taken from the coordinate input file if
    [`irest`][configs.workflow.amber.inputs.AmberInputsBase.irest] is `1`.
    """

    nrespa: int = Field(default=1, ge=1)
    r"""
    This variable allows the user to evaluate slowly-varying terms in the force field
    less frequently. For PME, "slowly-varying" (now) means the reciprocal sum.
    For generalized Born runs, the "slowly-varying" forces are those involving
    derivatives with respect to the effective radii, and pair interactions whose
    distances are greater than the "inner" cutoff, currently hard-wired at 8 Å.

    If [`nrespa`][configs.workflow.amber.inputs.AmberInputsBase.nrespa] > `1`
    these slowly-varying forces are evaluated every
    [`nrespa`][configs.workflow.amber.inputs.AmberInputsBase.nrespa] steps.
    The forces are adjusted appropriately, leading to an impulse at that step.
    If [`nrespa`][configs.workflow.amber.inputs.AmberInputsBase.nrespa]
    $\times$ [`dt`][configs.workflow.amber.inputs.AmberInputsBase.dt]
    is less than or equal to `4` fs then the energy conservation is not seriously
    compromised. However if
    [`nrespa`][configs.workflow.amber.inputs.AmberInputsBase.nrespa] $\times$
    [`dt`][configs.workflow.amber.inputs.AmberInputsBase.dt] > `4` fs then the
    simulation becomes less stable. Note that energies and related quantities are
    only accessible every
    [`nrespa`][configs.workflow.amber.inputs.AmberInputsBase.nrespa] steps, since
    the values at other times are meaningless.
    """

    temp0les: int = Field(default=-1, ge=-1)
    """This is the target temperature for all LES particles. If
    [`temp0les`][configs.workflow.amber.inputs.AmberInputsBase.temp0les] < `0`,
    a single temperature bath is used for all atoms, otherwise separate thermostats
    are used for LES and non-LES particles. Default is -1, corresponding to a single
    (weak-coupling) temperature bath.
    """

    tautp: float = Field(default=1.0, gt=0)
    """
    Time constant, in ps, for heat bath coupling for the system, if
    [`ntt`][configs.workflow.amber.inputs.AmberInputsBase.ntt] is `1`.
    Generally, values for
    [`tautp`][configs.workflow.amber.inputs.AmberInputsBase.tautp] should be in the
    range of `0.5` to `5.0` ps, with a smaller value providing tighter coupling to the
    heat bath and, thus, faster heating and a less natural trajectory. Smaller values of
    [`tautp`][configs.workflow.amber.inputs.AmberInputsBase.tautp] result in
    smaller fluctuations in kinetic energy, but larger fluctuations in the total energy.
    Values much larger than the length of the simulation result in a return
    to constant energy conditions.
    """

    vrand: int = Field(default=1000, ge=1)
    """
    If [`vrand`][configs.workflow.amber.inputs.AmberInputsBase.vrand] > `0`
    and [`ntt`][configs.workflow.amber.inputs.AmberInputsBase.ntt] is `2`, the
    velocities will be randomized to temperature
    [`temp0`][configs.workflow.amber.inputs.AmberInputsBase.temp0] every
    [`vrand`][configs.workflow.amber.inputs.AmberInputsBase.vrand] steps.
    """

    vlimit: float = Field(default=20.0, ge=0.0)
    """
    If not equal to `0.0`, then any component of the velocity that is greater than
    abs(VLIMIT) will be reduced to VLIMIT (preserving the sign). This can be used to
    avoid occasional instabilities in molecular dynamics runs. VLIMIT should
    generally be set to a value like 20 (the default), which is well above the most
    probable velocity in a Maxwell-Boltzmann distribution at room temperature.
    A warning message will be printed whenever the velocities are modified.
    Runs that have more than a few such warnings should be carefully examined.
    """

    nkija: int = Field(default=1, ge=1)
    r"""
    For use with [`ntt`][configs.workflow.amber.inputs.AmberInputsBase.ntt] if `9` or 10.
    For [`ntt`][configs.workflow.amber.inputs.AmberInputsBase.ntt] is `9`, this the
    number of substeps of [`dt`][configs.workflow.amber.inputs.AmberInputsBase.dt]
    when integrating the thermostat equations of motion, for greater accuracy.
    For [`ntt`][configs.workflow.amber.inputs.AmberInputsBase.ntt] is `10`, this
    specifies the number of additional auxiliary velocity variables v1 and v2, which
    will total `nkija` $\times$ v1 + `nkija` $\times$ v2
    """

    sinrtau: float = Field(default=1.0, gt=0.0)
    """
    For the SINR (Stochastic Isokinetic Nose-Hoover RESPA) integrator
    ([`ntt`][configs.workflow.amber.inputs.AmberInputsBase.ntt] is `10`), this
    specifies the time scale for determining the masses associated with the two
    auxiliary velocity variables v1 and v2 (e.g. thermostat velocities) and hence
    the magnitude of the coupling of the physical velocities with the auxiliary
    velocities. Generally this should be related to the time scale of the system.
    """

    baroscalingdir: Literal[0, 1, 2, 3] = Field(default=0)
    """
    Flag for pressure scaling direction control. Applicable when using Monte Carlo
    barostat ([`barostat`][configs.workflow.amber.inputs.AmberInputsBase.barostat]
    is `2`) with anisotropic pressure scaling
    ([`ntp`][configs.workflow.amber.inputs.AmberInputsBase.ntp] is `2`).

    <hr>

    **`0`**

    box size scales randomly ($x$, $y$ or $z$) each scaling step

    <hr>

    **`1`**

    box scales only along $x$-direction, dimensions along $y$-, $z$-axes are fixed

    <hr>

    **`2`**

    box scales only along $y$-direction, dimensions along $x$-, $z$-axes are fixed

    <hr>

    **`3`**

    box scales only along $z$-direction, dimensions along $x$-, $y$-axes are fixed
    """

    csurften: Literal[0, 1, 2, 3] = Field(default=0)
    """
    Flag for constant surface tension dynamics.

    <hr>

    **`0`**

    No constant surface tension.

    <hr>

    **`1`**

    Constant surface tension with interfaces in the $yz$ plane

    <hr>

    **`2`**

    Constant surface tension with interfaces in the $xz$ plane

    <hr>

    **`3`**

    Constant surface tension with interfaces in the $xy$ plane
    """

    gamma_ten: float = Field(default=0.0, ge=0.0)
    """Surface tension value in units of dyne/cm."""

    ninterface: int = Field(default=2, ge=2, le=9)
    """Number of interfaces in the periodic box. There must be at least two
    interfaces in the periodic box. Two interfaces is appropriate for a lipid
    bilayer system and is the default value.
    """

    tol: float = Field(default=0.00001, gt=0.0)
    """Relative geometrical tolerance for coordinate resetting in shake.
    Recommended maximum: < `0.00005` Angstroms.
    """

    jfastw: Literal[0, 4] = Field(default=0)
    """Fast water definition flag. By default, the system is searched for water
    residues, and special routines are used to SHAKE these systems.

    <hr>

    **`0`**

    Waters are identified by the default names, unless they are
    redefined, as described below.

    <hr>

    **`4`**

    Do not use the fast SHAKE routines for waters.

    !!! info
        The following variables allow redefinition of the default residue and atom names used by the
        program to determine which residues are waters.

        -   `WATNAM` The residue name the program expects for water. Default `WAT`.
        -   `OWTNM` The atom name the program expects for the oxygen of water.
            Default `O`.
        -   `HWTNM1` The atom name the program expects for the 1st H of water.
            Default `H1`.
        -   `HWTNM2` The atom name the program expects for the 2nd H of water.
            Default `H2`.
    """

    noshakemask: str = Field(default="")
    """
    String that specifies atoms that are not to be shaken
    (assuming that [`ntc`][configs.workflow.amber.inputs.AmberInputsBase.ntc]>1).
    Any bond that would otherwise be shaken by virtue of the
    [`ntc`][configs.workflow.amber.inputs.AmberInputsBase.ntc] flag, but which
    involves an atom flagged here, will **not** be shaken.
    Default is an empty string, which matches nothing. A typical use would be to remove
    SHAKE constraints from all or part of a solute, while still shaking rigid water
    models like TIPnP or SPC/E. Another use would be to turn off SHAKE constraints for
    the parts of the system that are being changed with thermodynamic integration, or
    which are the EVB or quantum regions of the system.

    If this option is invoked, then all parts of the potential must be evaluated,
    that is, [`ntf`][configs.workflow.amber.inputs.AmberInputsBase.ntf] must be `1`.
    The code enforces this by setting
    [`ntf`][configs.workflow.amber.inputs.AmberInputsBase.ntf] to `1` when a
    noshakemask string is present in the input.

    If you want the noshakemask to apply to all or part of the water molecules,
    you must also set
    [`jfastw`][configs.workflow.amber.inputs.AmberInputsBase.jfastw] to `4`, to turn
    off the special code for water SHAKE. (If you are not shaking waters, you
    presumably also want to issue the "set default FlexibleWater on" command in LEaP;
    see that chapter for more information.)
    """
