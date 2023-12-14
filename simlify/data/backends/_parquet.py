from typing import Any

from collections.abc import Iterable

import pyarrow as pa
import pyarrow.parquet as pq


def tabular_parquet(
    path: str,
    data: dict[str, Any],
    schema_fields: Iterable[Iterable[Any]],
    **kwargs: dict[str, Any],
) -> pa.Table:
    """Initialize `pyarrow.table` and save parquet file."""
    tab = pa.table(data, schema=pa.schema(schema_fields), **kwargs)
    pq.write_table(table=tab, where=path, **kwargs)
    return tab
