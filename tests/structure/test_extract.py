import os

from conftest import download_pdb

from simlify.structure import extract_atoms


def cache_1haj(dir_test: str) -> str:
    pdb_id = "1haj"
    path_pdb = os.path.join(dir_test, f"tmp/{pdb_id}.pdb")
    download_pdb(pdb_id, path_save=path_pdb)
    return path_pdb


def test_extract_everything(dir_test):
    """Extract everything out of NMR structure"""
    path_pdb = cache_1haj(dir_test)
    atoms = extract_atoms(
        path_topo=path_pdb, kwargs_universe={"topology_format": "pdb"}
    )

    assert atoms.n_atoms == 1310
    assert len(atoms.universe.trajectory) == 10


def test_extract_atoms(dir_test):
    """Extract chain A out of NMR structure"""
    path_pdb = cache_1haj(dir_test)
    atoms = extract_atoms(
        path_topo=path_pdb,
        select="chainID A",
        kwargs_universe={"topology_format": "pdb"},
    )

    assert atoms.n_atoms == 1085
    assert len(atoms.universe.trajectory) == 10
