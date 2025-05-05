# scaffold.py

import os

folders = [
    "scan_scn",
    "bulk/config",
    "bulk/scripts",
    "singlecell/config",
    "singlecell/scripts",
    "scripts",  # general-purpose (e.g. clean_inputs.py)
    "shared/signature_sets"
]

files = {
    "scan_scn/__init__.py": "",

    "bulk/Snakefile": "# Snakemake pipeline for bulk data\n",
    "bulk/config/config.yaml": "# Configuration for bulk pipeline\n",
    "bulk/scripts/score_scn.R": "# SCN scoring logic in R for bulk RNA-seq\n",

    "singlecell/Snakefile": "# Snakemake pipeline for single-cell data\n",
    "singlecell/config/config.yaml": "# Configuration for single-cell pipeline\n",
    "singlecell/scripts/score_scn.R": "# SCN scoring logic in R for scRNA-seq\n",

    "scripts/clean_inputs.py": "# Optional: input wrangling in Python\n",

    "shared/signature_sets/scn_markers.tsv": "geneA\ngeneB\n",

    ".gitignore": "data/\nresults/\n*.log\n__pycache__/\n"
}

for folder in folders:
    os.makedirs(folder, exist_ok=True)

for path, content in files.items():
    # Skip README.md and LICENSE if they already exist
    if os.path.exists(path):
        continue
    with open(path, "w") as f:
        f.write(content)

print("Project scaffold created.")