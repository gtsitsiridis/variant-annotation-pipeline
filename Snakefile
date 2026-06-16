"""variant_annotation — dataset-agnostic functional annotation of a variant set.

Tier 1 (VEP + plugins) + Tier 2 (AbSplice), assembled into output/annotations.parquet.
Configure paths in config.yaml (or --configfile config.local.yaml).
"""
from pathlib import Path

configfile: "config.yaml"

OUT = Path(config["output_dir"])
CHROMS = config["chromosomes"]
ABSPLICE = config.get("absplice", {}).get("enabled", False)

wildcard_constraints:
    chrom = r"chr[0-9XY]+",

include: "workflow/rules/vep.smk"
include: "workflow/rules/absplice.smk"
include: "workflow/rules/merge.smk"


rule all:
    input:
        OUT / "annotations.parquet",
