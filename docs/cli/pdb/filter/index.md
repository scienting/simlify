# `filter`

## Arguments

### `--records`

Every PDB file contains [many records](https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html) for depositors to use.
The main ones we care about are the ones in the [coordinate section](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html) including: [`ATOM`](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM), [`HETATM`](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#HETATM), [`TER`](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#TER), [`MODEL`](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#MODEL), [`ENDMDL`](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ENDMDL)
We can specify the records we want to keep using the `--records` argument.

```bash
--records ATOM,HETATM,TER,MODEL,ENDMDL 
```

!!! warning

    Records should always be separated by commas and there should no spaces between the records.
    For example, `--records ATOM, HETATM` and `--records ATOM;HETATM` will fail.

By default, not specifying this option will result in `ATOM`, `HETATM`, `TER`, `MODEL`, and `ENDMDL` being kept. 

### `--output`

Specifies where to write the new PDB file with only the specified records.

### `pdb_path`

Path to the PDB file to filter.

## Examples

Here are some examples of filtering lines from PDB files.

### Keeping coordinate information

```bash
simlify pdb filter --output filtered.pdb --records ATOM HETATM TER MODEL original.pdb
```
