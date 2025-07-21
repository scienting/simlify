# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

import os

from conftest import TMP_DIR

from simlify.schemas.amber import Amber22CLI, Amber22Schema


def test_render_amber_to_yaml():
    schema_amber = Amber22Schema()
    yaml_path = os.path.join(TMP_DIR, "amber22.yml")
    schema_amber.to_yaml(yaml_path)


def test_render_amber_write_input():
    schema_amber = Amber22Schema()
    input_path = os.path.join(TMP_DIR, "amber22.in")
    schema_amber.inputs.write_render(input_path)


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
