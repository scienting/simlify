from atomea.schemas import YamlIO
from pydantic import BaseModel


class ForcefieldSchemaBase(BaseModel, YamlIO):
    dna: str | None = None
    """Molecular mechanics force fields for DNA."""

    glycam: str | None = None
    """Molecular mechanics force fields for sugars"""

    ions: str | None = None
    """Molecular mechanics force fields for ions."""

    lipid: str | None = None
    """Molecular mechanics force fields for lipids."""

    protein: str | None = None
    """Molecular mechanics force field used to describe polypeptides."""

    rna: str | None = None
    """Molecular mechanics force fields for RNA."""

    small_molecule: str | None = None
    """Molecular mechanics force fields for small molecules."""

    water: str | None = None
    """Molecular mechanics force fields for water."""
