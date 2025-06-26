# PDB cleaning

## Trimming

???+ note "Prelude"

    ```text
    --8<-- "docs/case-studies/hiv-protease/files/structures/2PC0.pdb::53"
    ```

???+ note "Temperature factors"

    ```text
    --8<-- "docs/case-studies/hiv-protease/files/structures/2PC0.pdb:55:58"
    ```

## Removing water molecules

```text
--8<-- "docs/case-studies/hiv-protease/files/structures/2PC0.pdb:1735:1738"
```

```bash
--8<-- "docs/case-studies/hiv-protease/scripts/01-clean-pdb.sh:27:27"
```

## Removing non-protein molecules

```bash
--8<-- "docs/case-studies/hiv-protease/scripts/01-clean-pdb.sh:28:29"
```

## Removing duplicate atom lines

```bash
--8<-- "docs/case-studies/hiv-protease/scripts/01-clean-pdb.sh:32:52"
```

## Turning models into chains

```bash
--8<-- "docs/case-studies/hiv-protease/scripts/01-clean-pdb.sh:55:56"
```

## Renumber

```bash
--8<-- "docs/case-studies/hiv-protease/scripts/01-clean-pdb.sh:58:58"
```

## Result

<div id="hiv-protease-clean-view" class="mol-container"></div>

<script>
    document.addEventListener('DOMContentLoaded', (event) => {
        const viewer = molstar.Viewer.create('hiv-protease-clean-view', {
            layoutIsExpanded: false,
            layoutShowControls: false,
            layoutShowRemoteState: false,
            layoutShowSequence: true,
            layoutShowLog: false,
            layoutShowLeftPanel: false,
            viewportShowExpand: true,
            viewportShowSelectionMode: true,
            viewportShowAnimation: false,
            pdbProvider: 'rcsb',
        }).then(viewer => {
            // viewer.loadStructureFromUrl('./files/structures/2PC0-cleaned.pdb', format='pdb');
            viewer.loadSnapshotFromUrl("/case-studies/hiv-protease/files/molstar/2pc0-stripped.molx", "molx");
        });
    });
</script>


