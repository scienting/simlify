from dataclasses import dataclass

from pint import Quantity

from .. import ureg


@dataclass
class Coordinates:
    """Cartesian coordinates of a physical system."""

    dtype: str = "f8"

    ctype: str = "array"

    _units: Quantity = 1.0 * ureg.angstrom

    @property
    def units(self):
        """Default: Angstrom"""
        return self._units

    @units.setter
    def units(self, new_unit: Quantity) -> None:
        self._units *= new_unit
