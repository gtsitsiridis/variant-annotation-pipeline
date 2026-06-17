"""variant_annotation — dataset-agnostic functional annotation of a variant set.

VEP + plugins, NMD-Scanner, and AbSplice, assembled into output/annotations.parquet.
Configure paths in config.yaml (or --configfile config.local.yaml).
"""
from pathlib import Path

configfile: "config.yaml"

OUT = Path(config["output_dir"])
ABSPLICE = config.get("absplice", {}).get("enabled", False)
NMD = config.get("nmd", {}).get("enabled", False)

wildcard_constraints:
    chunk = r"\d+",

include: "workflow/rules/vep.smk"
include: "workflow/rules/nmd.smk"
include: "workflow/rules/absplice.smk"
include: "workflow/rules/merge.smk"


rule all:
    input:
        OUT / "annotations.parquet",
