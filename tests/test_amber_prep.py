import tempfile

from simlify.simulation.amber.topo import AmberTopoGen
from simlify.simulation.topo import run_gen_topo


class TestPrelimInfo:
    def test_1jc0(self, path_1jc0_prepped, amber_protein_standard_context):
        tleap_info = AmberTopoGen.dry_run(
            path_1jc0_prepped, amber_protein_standard_context, dir_work=None
        )
        assert len(tleap_info["duplicate_atoms"]) == 13
        assert tleap_info["unknown_residues"][-1]["residue_name"] == "CRO"
        assert tleap_info["n_water_molecules"] == 8761
        assert tleap_info["system_charge"] == -6.0

    def test_1jc0_with_cro(
        self,
        path_1jc0_prepped,
        amber_protein_standard_context,
        path_cro_fcrmod,
        path_cro_lib,
    ):
        amber_protein_standard_context.extra_lines_topo_gen = [
            'addAtomTypes { {"cc" "C" "sp2"} {"cd" "C" "sp2"} {"cf" "C" "sp2"} '
            '{"c" "C" "sp2"} {"nd" "N" "sp2"} {"nc" "N" "sp2"}{"ne" "N" "sp2"}'
            '{"nf" "N" "sp2"}{"ha" "H" "sp3"}{"oh" "O" "sp3"} }',
            f"xFPparams = loadamberparams {path_cro_fcrmod}",
            f"loadOff {path_cro_lib}",
        ]
        tleap_info = AmberTopoGen.dry_run(
            path_1jc0_prepped, amber_protein_standard_context
        )
        assert len(tleap_info["duplicate_atoms"]) == 13
        assert len(tleap_info["unknown_residues"]) == 0
        assert tleap_info["n_water_molecules"] == 8761
        assert tleap_info["system_charge"] == -7.0


class TestPrep:
    def test_1jc0_with_cro(
        self,
        path_1jc0_prepped,
        amber_protein_standard_context,
        path_cro_fcrmod,
        path_cro_lib,
    ):
        amber_protein_standard_context.extra_lines_topo_gen = [
            'addAtomTypes { {"cc" "C" "sp2"} {"cd" "C" "sp2"} {"cf" "C" "sp2"} '
            '{"c" "C" "sp2"} {"nd" "N" "sp2"} {"nc" "N" "sp2"}{"ne" "N" "sp2"}'
            '{"nf" "N" "sp2"}{"ha" "H" "sp3"}{"oh" "O" "sp3"} }',
            f"xFPparams = loadamberparams {path_cro_fcrmod}",
            f"loadOff {path_cro_lib}",
        ]
        path_topo = tempfile.NamedTemporaryFile().name
        path_coord = tempfile.NamedTemporaryFile().name
        run_gen_topo(
            path_1jc0_prepped,
            path_topo,
            path_coord,
            "simlify.simulation.amber.topo.AmberTopoGen",
            amber_protein_standard_context,
            dir_work=None,
        )
        with open(path_topo, "r", encoding="utf-8") as f:
            for line in f:
                if r"%FORMAT(10I8) " in line:
                    line = next(f)
                    assert line.split()[0] == "36341"
                    line = next(f)
                    assert line.split()[-2] == "46"
                    line = next(f)
                    assert line.split()[-3] == "1"
                    break
        with open(path_coord, "r", encoding="utf-8") as f:
            for line in f:
                if "default_name" in line:
                    line = next(f)
                    assert line.split()[0] == "36341"
                    line = next(f)
                    assert line.split()[-2] == "26.6466682"
                    break
