# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.


import numpy as np

from simlify.structure.dims import get_box_lengths, get_box_vectors, get_box_volume


def test_1haj_box_lengths(u_1haj):
    """Test that box lengths are computed correctly for 1HAJ."""
    pos = u_1haj.atoms.positions
    expected = np.array([43.5780003, 36.591, 27.844])
    np.testing.assert_allclose(get_box_lengths(pos), expected)


def test_1haj_box_vectors(u_1haj):
    """Test that box vectors are computed correctly for 1HAJ."""
    pos = u_1haj.atoms.positions
    expected = np.array(
        [[43.57800293, 0.0, 0.0], [0.0, 36.5909996, 0.0], [0.0, 0.0, 27.84399986]]
    )
    np.testing.assert_allclose(get_box_vectors(pos), expected)


def test_1haj_box_volume(u_1haj):
    """Test that box volume is computed correctly for 1HAJ."""
    pos = u_1haj.atoms.positions
    expected = 44399.00326322956
    assert np.isclose(get_box_volume(pos), expected)
