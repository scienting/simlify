from typing import Any

import numpy.typing as npt
import zarr


def array_zarr(
    path: str,
    arr: npt.NDArray[Any],
    mode: str = "a",
    dtype: str = "f8",
    **kwargs: dict[str, Any],
) -> zarr.core.Array:
    """Thin wrapper around `zarr.open`."""
    z1 = zarr.open(path, mode=mode, shape=arr.shape, dtype=dtype, **kwargs)
    z1[:] = arr
    return z1
