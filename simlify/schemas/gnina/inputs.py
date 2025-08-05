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
    flex: bool | None = Field(default=None)
    """Flag for using flexible docking.

    <hr>

    **`True`**

    <hr>
    
    Use gnina's flexible docking on flexible side chains, if any (PDBQT).

    <hr>

    **`False`**

    <hr>

    Do not use gnina's flexible docking procedure.
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
    <hr>
    `ad4_scoring`
    <hr>
    <hr>
    `default`
    <hr>
    <hr>
    `dkoes_fast`
    <hr>
    <hr> 
    `dkoes_scoring`
    <hr>
    <hr>
    `dkoes_scoring_old`
    <hr>
    <hr>
    `vina`
    <hr>
    <hr>
    `vinardo`
    <hr>
    """

    customr_scoring: str | None = Field(default=None)
    """Custom scoring function file.
    """

    custom_atoms: str | None = Field(default=None)
    """Custom atom type parameter file.
    """

    score_only: str | None = Field(default=None)
    """Score provided ligand pose
    """

    """
    The following control Gnina's Convolutional Neural Net (CNN) scoring.
    """

    """
    The following control Gnina's optional miscellaeneous parameters.
    """
