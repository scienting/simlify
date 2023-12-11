"""Generate the code reference pages and navigation."""

import os
from pathlib import Path

import mkdocs_gen_files

src = Path(__file__).parent.parent / "simulify"
write_dir = "api"

for path in sorted(src.rglob("*.py")):
    module_path = path.relative_to(src).with_suffix("")
    doc_path = path.relative_to(src).with_suffix(".md")
    full_doc_path = os.path.join(write_dir, doc_path)
    parts = tuple(module_path.parts)

    if parts[-1] == "__init__":
        parts = parts[:-1]
        full_doc_path = full_doc_path.replace("__init__.md", "index.md")
    elif parts[-1] == "__main__":
        continue

    if len(parts) == 0:
        parts = ("simulify",)

    with mkdocs_gen_files.open(full_doc_path, "w") as fd:
        identifier = ".".join(parts)
        print("::: " + identifier, file=fd)

    mkdocs_gen_files.set_edit_path(full_doc_path, path)
