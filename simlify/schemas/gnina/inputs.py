from atomea.schemas.io import YamlIO
from atomea.schemas.render import Render
from loguru import logger
from pydantic import BaseModel, Field


class GninaInputsBase(BaseModel, YamlIO, Render):
    """
    Gnina Input Base class.
    """

    """
    The following parameters below control Gnina's flexible docking schema.    
    """
    flex: str | None = Field(default=None)
    """Flag for using flexible docking.
    Use gnina's flexible docking on flexible side chains, if any (PDBQT).
    """

    flexres: list[str] | None = Field(default=None)
    """User parameter to set flexible sidechains as a list of chain residue
    key value pair in comma seperated format.
    Flexible side chains specified by comma separated list of `chain:resid`.
    """

    flexdist_ligand: str | None = Field(default=None)
    """User specified ligand file for using flexible sidechains.
    Ligand to use for `flexdist`.
    """

    flexdist: float | None = Field(default=None)
    """User specified distance to use from the ligand in `flexdist_ligand`.
    Sets all side chains within specified distance to `flexdist_ligand` to flexible.
    """

    flex_limit: int | None = Field(default=None)
    """Hard limit for the number of flexible residues.
    """

    flex_max: int | None = Field(default=None)
    """Retain at at most the closest `flex_max` flexible residues.
    """

    """
    The following parameters control Gnina's serach space.
    """

    center_x: float | None = Field(default=None)
    """Defines the center of the docking box.
    X coordinate of the center.
    """

    center_y: float | None = Field(default=None)
    """Defines the center of the docking box.
    Y coordinate of the center.
    """

    center_z: float | None = Field(default=None)
    """Defines the center of the docking box.
    Z coordinate of the center.
    """

    size_x: float | None = Field(default=None)
    """Defines the size of the docking box.
    Size in the X dimension (Å).
    """

    size_y: float | None = Field(default=None)
    """Defines the size of the docking box.
    Size in the Y dimension (Å).
    """

    size_z: float | None = Field(default=None)
    """Defines the size of the docking box.
    Size in the Z dimension (Å).
    """

    autobox_ligand: str | None = Field(default=None)
    """Name of the ligand file to define the autobox around.
    Ligand to use for autobox. A multi-ligand file still only defines a single box.
    """

    autobox_add: float = Field(default=4.0)
    """Amount of buffer space to add to auto-generated box
    (default +4 on all six sides)    
    """

    autobox_extend: float = Field(default=1)
    """Expand the autobox if needed to ensure the input conformation of the ligand
    being docked can freely rotate within the box.    
    """

    no_lig: bool | None = Field(default=None)
    """Specify no ligand; for sampling/minimizing flexible residues
    """

    """
    The following control Gnina's covalent docking parameters.
    """

    covalent_rec_atom: str | list[float] | None = Field(default=None)
    """Receptor atom ligand is covalently bound to. Can be specified as 
    chain:resnum:atom_name or as x,y,z Cartesian coordinates.
    """

    covalent_lig_atom_pattern: str | None = Field(default=None)
    """SMARTS expression for ligand atom that will covalently bind protein.        
    """

    covalent_lig_atom_position: list[float] | None = Field(default=None)
    """Optional: Initial placement of covalently bonding ligand atom in x,y,z Cartesian 
    coordinates.  If not specified, OpenBabel's GetNewBondVector function will
    be used to position ligand.
    """

    covalent_fix_lig_atom_position: bool | None = Field(default=None)
    """If `covalent_lig_atom_position` is specified, fix the ligand atom to this 
    position as opposed to using this position to define the initial structure.
    """

    covalent_bond_order: int | None = Field(default=None)
    """Bond order of covalent bond. Default 1.
    """

    covalent_optimize_lig: bool | None = Field(default=None)
    """Optimize the covalent complex of ligand and residue using UFF. This will change 
    bond angles and lengths of the ligand.    
    """

    """
    The following control Gnina's scoring and minimization parameters.
    """

    scoring: str | None = Field(default=None)
    """Specify alternative built-in scoring function:
    `ad4_scoring`
    `default`
    `dkoes_fast`
    `dkoes_scoring`
    `dkoes_scoring_old`
    `vina`
    `vinardo`
    """

    custom_scoring: str | None = Field(default=None)
    """Custom scoring function file.
    """

    custom_atoms: str | None = Field(default=None)
    """Custom atom type parameter file.
    """

    score_only: str | None = Field(default=None)
    """Score provided ligand pose
    """

    local_only: bool | None = Field(default=None)
    """Local search only using autobox (you probably want to use `minimize`)
    """

    minimize: bool | None = Field(default=None)
    """Energy minimization.       
    """

    randomize_only: bool | None = Field(default=None)
    """Generate random poses, attempting to avoid clashes. 
    """

    num_mc_steps: int | None = Field(default=None)
    """Fixed number of monte carlo steps to take in each chain.
    """

    max_mc_steps: int | None = Field(default=None)
    """Cap on number of monte carlo steps to take in each chain.
    """

    num_mc_saved: int | None = Field(default=None)
    """Number of top poses saved in each monte carlo chain.
    """

    temperature: float | None = Field(default=None)
    """Temperature for metropolis accept criterion.
    """

    minimize_iters: float | None = Field(default=None)
    """number iterations of steepest descent; default scales with
    rotors and usually isn't sufficient for convergence.    
    """

    accurate_line: bool | None = Field(default=None)
    """Use accurate line search.
    """

    simple_ascent: bool | None = Field(default=None)
    """Use simple gradient ascent.
    """

    minimize_early_term: bool | None = Field(default=None)
    """Stop minimization before convergence conditions are fully met.
    """

    minimize_single_full: bool | None = Field(default=None)
    """During docking perform a single full minimization instead of
    a truncated pre-evaluate followed by a full.
    """

    approximation: str | None = Field(default=None)
    """Approximation to use:
    `linear`
    `spline`
    `exact`
    """

    factor: float | None = Field(default=None)
    """Approximation factor: higher results in a finer-grained approximation.
    """

    force_cap: float | None = Field(default=None)
    """Max allowed force; lower values more gently minimize clashing structures. 
    """

    user_grid: str | None = Field(default=None)
    """Autodock map file for user grid data based calculations.
    """

    user_grid_lambda: float | None = Field(default=None)
    """Scales user_grid and functional scoring.        
    """

    print_terms: bool | None = Field(default=None)
    """Print all available terms with default parameterizations.
    """

    print_atom_types: bool | None = Field(default=None)
    """Print all available atom types.
    """

    """
    The following control Gnina's Convolutional Neural Net (CNN) scoring.
    """

    cnn_scoring: str | None = Field(default=None)
    """Amount of CNN scoring:
    `none`
    `rescore`
    `(default)`
    `refinement`
    `metrorescore`
    `metropolis+rescore`
    `metrorefine` 
    `metropolis+refine`
    `all`
    """

    cnn: str | None = Field(default=None)
    """Built-in model to use, specify PREFIX_ensemble to evaluate an ensemble of
    models starting with PREFIX:
    `all_default_to_default_1_3_1`
    `all_default_to_default_1_3_2` 
    `all_default_to_default_1_3_3` 
    `crossdock_default2018` 
    `crossdock_default2018_1` 
    `crossdock_default2018_1_3` 
    `crossdock_default2018_1_3_1` 
    `crossdock_default2018_1_3_2` 
    `crossdock_default2018_1_3_3` 
    `crossdock_default2018_1_3_4` 
    `crossdock_default2018_2` 
    `crossdock_default2018_3` 
    `crossdock_default2018_4` 
    `crossdock_default2018_KD_1` 
    `crossdock_default2018_KD_2` 
    `crossdock_default2018_KD_3` 
    `crossdock_default2018_KD_4` 
    `crossdock_default2018_KD_5 default1.0` 
    `default2017 dense dense_1 dense_1_3` 
    `dense_1_3_1 dense_1_3_2 dense_1_3_3` 
    `dense_1_3_4 dense_1_3_PT_KD` 
    `dense_1_3_PT_KD_1 dense_1_3_PT_KD_2` 
    `dense_1_3_PT_KD_3 dense_1_3_PT_KD_4` 
    `dense_1_3_PT_KD_def2018` 
    `dense_1_3_PT_KD_def2018_1` 
    `dense_1_3_PT_KD_def2018_2` 
    `dense_1_3_PT_KD_def2018_3` 
    `dense_1_3_PT_KD_def2018_4 dense_2 dense_3` 
    `dense_4 fast general_default2018` 
    `general_default2018_1` 
    `general_default2018_2` 
    `general_default2018_3` 
    `general_default2018_4` 
    `general_default2018_KD_1` 
    `general_default2018_KD_2` 
    `general_default2018_KD_3` 
    `general_default2018_KD_4` 
    `general_default2018_KD_5` 
    `redock_default2018 redock_default2018_1` 
    `redock_default2018_1_3` 
    `redock_default2018_1_3_1` 
    `redock_default2018_1_3_2` 
    `redock_default2018_1_3_3` 
    `redock_default2018_1_3_4` 
    `redock_default2018_2 redock_default2018_3` 
    `redock_default2018_4 redock_default2018_KD`
    ``_1 redock_default2018_KD_2` 
    `redock_default2018_KD_3` 
    `redock_default2018_KD_4` 
    `redock_default2018_KD_5`        
    """

    cnn_model: str | None = Field(default=None)
    """Torch cnn model file; if not specified a default model ensemble will be used.
    """

    cnn_rotation: int | None = Field(default=None)
    """Evaluate multiple rotations of pose (max `24`).
    """

    cnn_mix_emp_force: bool | None = Field(default=None)
    """Merge CNN and empirical minus forces.
    """

    cnn_mix_emp_energy: bool | None = Field(default=None)
    """Merge CNN and empirical energy.
    """

    cnn_empirical_weight: float | None = Field(default=None)
    """Weight for scaling and merging empirical force and energy.
    """

    cnn_center_x: float | None = Field(default=None)
    """X coordinate of the CNN center.
    """

    cnn_center_y: float | None = Field(default=None)
    """Y coordinate of the CNN center.
    """

    cnn_center_z: float | None = Field(default=None)
    """Z coordinate of the CNN center.
    """

    cnn_verbose: bool | None = Field(default=None)
    """Enable verbose output for CNN debugging.
    """

    """
    The following control Gnina's optional miscellaeneous parameters.
    """

    out_flex: str | None = Field(default=None)
    """Output file for flexible receptor residues.
    """

    atom_terms: str | None = Field(default=None)
    """Optionally write per-atom interaction term values
    """

    atom_term_data: bool | None = Field(default=None)
    """Embedded per-atom interaction terms in output sd data.
    """

    pose_sort_order: str | None = Field(default=None)
    """How to sort docking results:
    `CNNscore` 
    `default`
    `CNNaffinity`
    `Energy`    
    """

    full_flex_output: bool | None = Field(default=None)
    """Output entire structure for out_flex, not just flexible residues.
    """

    cpu: int = Field(default=1)
    """The number of CPUs to use (the default is to try to detect the
    number of CPUs or, failing that, use 1)    
    """

    seed: int | None = Field(default=None)
    """Explicit random seed.
    """

    exhaustiveness: int = Field(default=8)
    """Exhaustiveness of the global search (roughly proportional to time).
    Defaults to `8`.
    """

    num_modes: int = Field(default=9)
    """Maximum number of binding modes to generate.
    """

    min_rmsd_filter: float = Field(default=1)
    """RMSD value used to filter final poses to remove redundancy.        
    """

    quiet: bool | None = Field(default=None)
    """Suppress output messages
    """

    addH: str | None = Field(default=None)
    """Automatically add hydrogens in ligands (on by default).
    """

    stripH: str | None = Field(default=None)
    """Remove polar hydrogens from molecule *after* performing
    atom typing for efficiency (off by default - nonpolar are always removed)
    """

    device: int | None = Field(default=0)
    """GPU device to use
    """

    no_gpu: bool | None = Field(default=None)
    """Disable GPU acceleration, even if available.
    """
