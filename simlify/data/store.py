import os

import MDAnalysis as mda
from loguru import logger

from .backends._parquet import tabular_parquet
from .backends._zarr import array_zarr


class SimStore:
    """Store simulation data into a single format."""

    def __init__(self, backend_array=array_zarr, backend_tabular=tabular_parquet):
        self.backend_array = backend_array
        self.backend_tabular = backend_tabular

    def process(self, dir_write, topology, *args, **kwargs):
        """Store trajectory according to backends."""
        u = mda.Universe(topology, *args, **kwargs)
        try:
            self.backend_array(
                os.path.join(dir_write, "coordinates"), u.atoms.positions
            )
        except mda.exceptions.NoDataError:
            logger.debug("No positions in universe")
