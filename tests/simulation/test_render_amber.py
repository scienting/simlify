import os

from conftest import TMP_DIR

from simlify.configs.amber import Amber22CLI, Amber22Config


def test_render_amber_to_yaml():
    config_amber = Amber22Config()
    yaml_path = os.path.join(TMP_DIR, "amber22.yml")
    config_amber.to_yaml(yaml_path)


def test_render_amber_write_input():
    config_amber = Amber22Config()
    input_path = os.path.join(TMP_DIR, "amber22.in")
    config_amber.inputs.write_render(input_path)


def test_render_amber_cli():
    amber_cli = Amber22CLI(
        mdin="mdin",
        mdout="mdout",
        mdinfo="mdinfo",
        prmtop="prmtop",
        inpcrd="inpcrd",
    )
    amber_command = amber_cli.render()

    assert len(amber_command) == 1
    assert amber_command[0] == "pmemd -i mdin -o mdout -inf mdinfo -p prmtop -c inpcrd"
