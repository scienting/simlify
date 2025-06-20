import os

from conftest import download_pdb

from simlify.structure.pdb.numbering.main import run_unify_numbering


def cache_5khb(dir_test: str) -> str:
    pdb_id = "5khb"
    path_pdb = os.path.join(dir_test, f"tmp/{pdb_id}.pdb")
    download_pdb(pdb_id, path_save=path_pdb)
    return path_pdb


def test_unify_5khb_noreset(dir_test):
    path_pdb = cache_5khb(dir_test)
    lines = run_unify_numbering(path_pdb, reset_initial_resid=False)
    line_hetatm = lines[140]
    line_hetatm_ref = (
        "HETATM    6  CG  FME A   1       0.765  -1.863  -1.806  1.00 35.55           C"
    )
    assert line_hetatm.strip() == line_hetatm_ref
    line_atom = lines[155]
    line_atom_ref = (
        "ATOM     21  CA  GLY A   2       2.192   3.731  -2.537  1.00 30.10           C"
    )
    assert line_atom.strip() == line_atom_ref


def test_unify_5khb_reset(dir_test):
    path_pdb = cache_5khb(dir_test)
    lines = run_unify_numbering(path_pdb, reset_initial_resid=True)
    line_hetatm = lines[140]
    line_hetatm_ref = (
        "HETATM    6  CG  FME A   1       0.765  -1.863  -1.806  1.00 35.55           C"
    )
    assert line_hetatm.strip() == line_hetatm_ref
    line_atom = lines[155]
    line_atom_ref = (
        "ATOM     21  CA  GLY A   2       2.192   3.731  -2.537  1.00 30.10           C"
    )
    assert line_atom.strip() == line_atom_ref


def test_multichain_toy_pdb(multichain_toy_pdb):
    """This tests the functionality of the function using a two chain toy pdb model with two identical chains"""
    path_pdb = str(multichain_toy_pdb)
    lines = run_unify_numbering(path_pdb, reset_initial_resid=False)
    line_first_chain_ref = "ATOM      1  N   MET A   1      11.104  13.207  10.123  1.00 20.00           N  \n"

    line_first_chain = lines[1]
    assert line_first_chain == line_first_chain_ref


def test_multichain_toy_pdb_reset(multichain_toy_pdb):
    """This tests the functionality of the function using a two chain toy pdb model with two identical chains and tests the reset resid boolean"""
    path_pdb = str(multichain_toy_pdb)
    lines = run_unify_numbering(path_pdb, reset_initial_resid=True)
    line_first_chain_ref = "ATOM      1  N   MET A   1      11.104  13.207  10.123  1.00 20.00           N  \n"

    line_first_chain = lines[1]
    assert line_first_chain == line_first_chain_ref


def test_multichain_hetatm_toy_pdb(multichain_hetatm_toy_pdb):
    """This tests the functionality of the function using a two chain with toy pdb model with two identical chains 'HETATM' records"""
    path_pdb = str(multichain_hetatm_toy_pdb)
    lines = run_unify_numbering(path_pdb, reset_initial_resid=False)
    line_first_chain_ref = "HETATM    9  N   ALA B   1      15.104  12.207  15.123  1.00 20.00           N  \n"
    line_first_chain = lines[10]
    assert line_first_chain == line_first_chain_ref


def test_multichain_hetatm_toy_pdb_reset(multichain_hetatm_toy_pdb):
    """This tests the functionality of the function using a two chain with toy pdb model with two identical chains 'HETATM' records and tests the boolean"""
    path_pdb = str(multichain_hetatm_toy_pdb)
    lines = run_unify_numbering(path_pdb, reset_initial_resid=True)
    line_first_chain_ref = "HETATM    9  N   ALA B   1      15.104  12.207  15.123  1.00 20.00           N  \n"

    line_first_chain = lines[10]
    assert line_first_chain == line_first_chain_ref


def test_seq_multichain_toy_pdb(multichain_seq_hetatm_toy_pdb):
    """This function tests the toy multichain model with sequential resids with 'HETATM' records"""
    path_pdb = str(multichain_seq_hetatm_toy_pdb)
    lines = run_unify_numbering(path_pdb, reset_initial_resid=False)
    line_first_chain_ref = "HETATM    9  N   ALA B   3      15.104  12.207  15.123  1.00 20.00           N  \n"
    line_first_chain = lines[10]
    for line in lines:
        print(line)
    assert line_first_chain == line_first_chain_ref


def test_seq_multichain_toy_pdb_reset(multichain_seq_hetatm_toy_pdb):
    """This function tests the boolean statement with a sequential resID multichain toy pdb model."""
    path_pdb = str(multichain_seq_hetatm_toy_pdb)
    lines = run_unify_numbering(path_pdb, reset_initial_resid=True)
    line_first_chain_ref = "HETATM    9  N   ALA B   1      15.104  12.207  15.123  1.00 20.00           N  \n"

    line_first_chain = lines[10]
    assert line_first_chain == line_first_chain_ref


def test_broken_multichain_toy_pdb(broken_multichain_toy_pdb):
    """This tests a completely broken pdb model missing termination records for chains and nonsequential resid numbering."""
    path_pdb = str(broken_multichain_toy_pdb)
    lines = run_unify_numbering(path_pdb, reset_initial_resid=False)
    line_first_chain_ref = "HETATM    9  N   ALA B   5      15.104  12.207  15.123  1.00 20.00           N  \n"
    line_first_chain = lines[10]
    for line in lines:
        print(line)
    assert line_first_chain == line_first_chain_ref


def test_broken_multichain_toy_pdb_reset(broken_multichain_toy_pdb):
    """This tests a completely broken pdb model missing termination records for chains and nonsequential resid number with resetting the resids"""
    path_pdb = str(broken_multichain_toy_pdb)
    lines = run_unify_numbering(path_pdb, reset_initial_resid=True)
    line_first_chain_ref = "HETATM    9  N   ALA B   1      15.104  12.207  15.123  1.00 20.00           N  \n"

    line_first_chain = lines[10]
    assert line_first_chain == line_first_chain_ref
