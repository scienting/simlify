# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

from simlify.simulation.amber.prep import AmberSimPrep


class TestAmberSimPrep:
    def test_amber_run_command_prep_01_min(self, amber_sim_standard_config):
        """Prepare run_command for amber simulation."""
        simlify_config = amber_sim_standard_config
        prep = AmberSimPrep()

        simlify_config = prep.prepare_sim_config(simlify_config)

        run_command = AmberSimPrep.get_stage_run_command(simlify_config)
        assert len(run_command) == 4
        assert run_command[0] == ""
        assert run_command[1] == "echo 'Starting 01_min'"
        assert run_command[2] == "date"

        amber_command = run_command[3]
        assert "None" not in amber_command
        assert "mpirun -np 8 pmemd.MPI -i" == amber_command[:25]
        amber_command_split = amber_command[24:].split(" -")
        assert "i " == amber_command_split[0][:2]
        assert "01_min.in" in amber_command_split[0]

        assert "o " == amber_command_split[1][:2]
        assert "01_min.out" in amber_command_split[1]

        assert "inf " == amber_command_split[2][:4]
        assert "01_min.mdinfo" in amber_command_split[2]

        assert "p " == amber_command_split[3][:2]
        assert "mol.prmtop" in amber_command_split[3]

        assert "c " == amber_command_split[4][:2]
        assert "mol.inpcrd" in amber_command_split[4]

        assert "x " == amber_command_split[5][:2]
        assert "01_min.nc" in amber_command_split[5]

        assert "r " == amber_command_split[6][:2]
        assert "01_min.rst" in amber_command_split[6]

    def test_amber_prepare_stage(self, amber_sim_standard_config):
        """Test preparing input file."""
        simlify_config = amber_sim_standard_config
        prep = AmberSimPrep()
        simlify_config = prep.prepare_sim_config(simlify_config)

        stage_input_lines, _ = AmberSimPrep.prepare_stage(simlify_config, write=False)
        stage_input_lines_ref = [
            "",
            "&cntrl",
            "    imin=0,",
            "    irest=1,",
            "    ntx=5,",
            "    ntmin=1,",
            "    maxcyc=9999,",
            "    ncyc=10,",
            "    ig=-1,",
            "    dt=0.002,",
            "    nstlim=500000,",
            "    nscm=500,",
            "    ntr=1,",
            "    restraint_wt=0.5,",
            "    restraintmask=\"!(:WAT) & (@C,CA,N,O,O5',P,O3',C3',C4',C5')\",",
            "    ntb=2,",
            "    ntf=2,",
            "    ntc=2,",
            "    cut=10.0,",
            "    ntt=3,",
            "    tempi=100.0,",
            "    temp0=300.0,",
            "    gamma_ln=5.0,",
            "    ntp=1,",
            "    barostat=2,",
            "    pres0=1.01325,",
            "    mcbarint=100,",
            "    comp=44.6,",
            "    taup=1.0,",
            "    ntxo=2,",
            "    ntwr=5000,",
            "    ntpr=500,",
            "    ntwx=5000,",
            "    ioutfm=1,",
            "    iwrap=1,",
            "    nmropt=0,",
            "    ntave=0,",
            "    ntwv=0,",
            "    ionstepvelocities=0,",
            "    ntwf=0,",
            "    ntwe=0,",
            "    ntwprt=0,",
            "    idecomp=0,",
            "    ibelly=0,",
            "    dx0=0.01,",
            "    drms=0.0001,",
            "    t=0.0,",
            "    nrespa=1,",
            "    temp0les=-1,",
            "    tautp=1.0,",
            "    vrand=1000,",
            "    vlimit=20.0,",
            "    nkija=1,",
            "    sinrtau=1.0,",
            "    baroscalingdir=0,",
            "    csurften=0,",
            "    gamma_ten=0.0,",
            "    ninterface=2,",
            "    tol=1e-05,",
            "    jfastw=0,",
            "&end",
            "",
        ]
        assert stage_input_lines == stage_input_lines_ref
