from simlify.structure.solvent import get_ion_counts


class TestNumberOfIons:
    def test_1jc0(
        self,
        amber_sim_standard_config,
    ):
        ion_counts = get_ion_counts(
            simlify_config=amber_sim_standard_config,
            charge_net=-7.0,
            n_waters=8761,
        )
        assert ion_counts["charge_anion_num"] == 23
        assert ion_counts["charge_cation_num"] == 30
