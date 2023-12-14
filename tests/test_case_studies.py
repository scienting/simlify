import os
import shutil

from simlify.simulation.contexts import SimContextManager
from simlify.simulation.prep import run_sim_slurm_prep
from simlify.simulation.topo import run_gen_topo


def test_amber_f7498a8c_d021_491c_a343_10151e81434a(dir_structures, dir_amber_sims):
    """Generate all simulation files for f7498a8c-d021-491c-a343-10151e81434a"""
    path_config = os.path.join(
        dir_amber_sims, "f7498a8c-d021-491c-a343-10151e81434a" + ".yml"
    )
    dir_write = os.path.join(dir_amber_sims, "f7498a8c-d021-491c-a343-10151e81434a")
    if os.path.exists(dir_write):
        shutil.rmtree(dir_write)
    os.makedirs(dir_write, exist_ok=False)
    topo_dir = os.path.join(dir_write, "topo")
    os.makedirs(topo_dir, exist_ok=False)

    sim_context_manager = SimContextManager(path_config)
    sim_context_manager.dir_write = dir_write

    # Topology
    path_structure = os.path.join(
        dir_structures, "f7498a8c-d021-491c-a343-10151e81434a.pdb"
    )
    path_topo_write = "topo/mol.prmtop"
    path_coord_write = "topo/mol.inpcrd"
    import_string = "simlify.simulation.amber.topo.AmberTopoGen"
    run_gen_topo(
        path_structure,
        path_topo_write,
        path_coord_write,
        import_string,
        sim_context_manager,
        dir_write,
    )
    # TODO: Add checks

    # Simulation files
    name_job = "f7498a8c-d021-491c-a343-10151e81434a"
    path_run_write = "run.sh"
    path_slurm_write = "submit.slurm"
    prep_class_string = "simlify.simulation.amber.prep.AmberSimPrep"

    run_sim_slurm_prep(
        name_job,
        dir_write,
        path_run_write,
        path_slurm_write,
        prep_class_string,
        sim_context_manager,
    )
