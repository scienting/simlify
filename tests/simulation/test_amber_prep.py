import os

from simlify.simulation.amber.topo import AmberTopoGen
from simlify.simulation.topo import run_gen_topo


class TestPrelimInfo:
    def test_1jc0(self, path_1jc0_prepped, amber_sim_standard_config):
        tleap_info = AmberTopoGen.dry_run(path_1jc0_prepped, amber_sim_standard_config)
        assert len(tleap_info["duplicate_atoms"]) == 0
        assert tleap_info["unknown_residues"][-1]["residue_name"] == "CRO"
        assert tleap_info["solvent_molecules_num"] == 8087
        assert tleap_info["charge_net"] == -8.0

    def test_1jc0_with_cro(
        self,
        path_1jc0_prepped,
        amber_sim_standard_config,
        path_cro_fcrmod,
        path_cro_lib,
    ):
        amber_sim_standard_config.topology.append_lines = [
            'addAtomTypes { {"cc" "C" "sp2"} {"cd" "C" "sp2"} {"cf" "C" "sp2"} '
            '{"c" "C" "sp2"} {"nd" "N" "sp2"} {"nc" "N" "sp2"} {"ne" "N" "sp2"}'
            '{"nf" "N" "sp2"} {"ha" "H" "sp3"} {"oh" "O" "sp3"} }',
            f"xFPparams = loadamberparams {path_cro_fcrmod}",
            f"loadOff {path_cro_lib}",
        ]
        tleap_info = AmberTopoGen.dry_run(path_1jc0_prepped, amber_sim_standard_config)
        assert len(tleap_info["duplicate_atoms"]) == 0
        assert len(tleap_info["unknown_residues"]) == 0
        assert tleap_info["solvent_molecules_num"] == 8087
        assert tleap_info["charge_net"] == -9.0


class TestPrep:
    def test_1jc0_with_cro(
        self,
        path_1jc0_prepped,
        amber_sim_standard_config,
        path_cro_fcrmod,
        path_cro_lib,
    ):
        amber_sim_standard_config.topology.append_lines = [
            'addAtomTypes { {"cc" "C" "sp2"} {"cd" "C" "sp2"} {"cf" "C" "sp2"} '
            '{"c" "C" "sp2"} {"nd" "N" "sp2"} {"nc" "N" "sp2"}{"ne" "N" "sp2"}'
            '{"nf" "N" "sp2"}{"ha" "H" "sp3"}{"oh" "O" "sp3"} }',
            f"xFPparams = loadamberparams {path_cro_fcrmod}",
            f"loadOff {path_cro_lib}",
        ]

        run_gen_topo(
            path_1jc0_prepped,
            "simlify.simulation.amber.topo.AmberTopoGen",
            amber_sim_standard_config,
        )

        path_topo = os.path.join(
            amber_sim_standard_config.run.dir_work,
            amber_sim_standard_config.run.dir_input,
            amber_sim_standard_config.engine.cli.prmtop,
        )
        with open(path_topo, "r", encoding="utf-8") as f:
            for line in f:
                if r"%FORMAT(10I8) " in line:
                    line = next(f)
                    assert line.split()[0] == "37381"
                    line = next(f)
                    assert line.split()[-2] == "46"
                    line = next(f)
                    assert line.split()[-3] == "1"
                    break

        path_coord = os.path.join(
            amber_sim_standard_config.run.dir_work,
            amber_sim_standard_config.run.dir_input,
            amber_sim_standard_config.engine.cli.inpcrd,
        )
        with open(path_coord, "r", encoding="utf-8") as f:
            for line in f:
                if "default_name" in line:
                    line = next(f)
                    assert line.split()[0] == "37381"
                    line = next(f)
                    assert line.split()[-2] == "23.6286682"
                    break
