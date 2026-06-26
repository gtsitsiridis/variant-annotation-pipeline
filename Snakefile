"""variant_annotation — dataset-agnostic functional annotation of a variant set.

A tiered funnel over a generic VCF (small variants ID=chrom_start_ref_alt; SVs ID=sv_id):

  Tier 0  fastVEP        broad per-transcript consequence/distance over the WHOLE set,
                         per track -> <track>/basic_annotations.parquet   (cheap, ALWAYS on)
             | select (deep.window) from the small track
  Deep    VEP + plugins  LOFTEE/CADD/SpliceAI/PrimateAI/AlphaMissense (+ NMD, AbSplice) on the
                         narrow analysis set (chunked) -> annotations.parquet
  E2G     ENCODE-rE2G    per-variant enhancer -> target gene + rE2G score -> e2g/e2g.parquet

Tier 0 always runs (only needs gff3+fasta). The deep VEP/NMD/AbSplice tier and E2G are opt-in
via config flags (deep.enabled / nmd.enabled / absplice.enabled / e2g.enabled), each needing
reference data. E2G is a STANDALONE parquet (grain variant x target_gene) joined by the
consumer on variant_id, NOT folded into the transcript-level annotations.parquet. Configure
paths in config.yaml (or --configfile config.local.yaml), or pass config as a Snakemake
`module`.
"""
from pathlib import Path

configfile: "config.yaml"

OUT = Path(config["output_dir"])
DEEP = config.get("deep", {}).get("enabled", False)
NMD = config.get("nmd", {}).get("enabled", False)
ABSPLICE = config.get("absplice", {}).get("enabled", False)
E2G = config.get("e2g", {}).get("enabled", False)

wildcard_constraints:
    chunk = r"\d+",

# Tier 0 — always (defines TRACKS, BASIC_SMALL).
include: "workflow/rules/fastvep.smk"

# Deep VEP/NMD/AbSplice tier — opt-in. select.smk defines ANALYSIS_SET (chunk_vcf reads it);
# merge spines on the chunked vep.parquet + canonical, folding in nmd/absplice when enabled.
if DEEP:
    include: "workflow/rules/select.smk"
    include: "workflow/rules/vep.smk"
    if NMD:
        include: "workflow/rules/nmd.smk"
    if ABSPLICE:
        include: "workflow/rules/absplice.smk"
    include: "workflow/rules/merge.smk"

if E2G:
    include: "workflow/rules/e2g.smk"


def _targets():
    t = [OUT / track / "basic_annotations.parquet" for track in TRACKS]   # Tier 0 always
    if DEEP:
        t.append(OUT / "annotations.parquet")           # deep merged (vep[+nmd+absplice])
    if E2G:
        t.append(OUT / "e2g" / "e2g.parquet")           # standalone (joined on variant_id)
    return [str(x) for x in t]


rule all:
    input:
        _targets(),
