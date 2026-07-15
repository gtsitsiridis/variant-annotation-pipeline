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

# ── Per-tool flags (the retired "deep tier" umbrella is gone) ─────────────────────────
# fastVEP is the always-on base. Every other tool has its own `enabled` + `distance`. Each tool
# writes OUT/<tool>/variant_type={SNV,indel,SV}/ keyed on variant_id. SV annotation is fastVEP-only
# for now (`fastvep.include_sv`); every other tool is small-variant only.
VEP = config.get("vep", {}).get("enabled", False)
NMD = config.get("nmd", {}).get("enabled", False)
ABSPLICE = config.get("absplice", {}).get("enabled", False)
E2G = config.get("e2g", {}).get("enabled", False)

# SV annotation is fastVEP-only for now. Reject include_sv on any other tool so the flag never lies.
for _t in ("vep", "nmd", "absplice", "e2g"):
    if config.get(_t, {}).get("include_sv", False):
        raise ValueError(f"{_t}.include_sv is not supported — SV annotation is fastVEP-only. Remove it.")

wildcard_constraints:
    chunk = r"\d+",
    track = r"small|sv",
    vtype = r"SNV|indel|SV",

include: "workflow/rules/fastvep.smk"        # Tier 0 base — always; defines FASTVEP_TARGETS

# NOTE: vep / nmd / absplice / e2g are being migrated to the per-tool <tool>/variant_type=<vtype>
# layout (carrying variant_id, own distance, include_sv). Until that lands, only fastVEP builds;
# enabling one of them fails loudly rather than running the old deep-tier machinery.
_PENDING = [t for t, on in (("vep", VEP), ("nmd", NMD), ("absplice", ABSPLICE), ("e2g", E2G)) if on]
if _PENDING:
    raise ValueError(
        "tools not yet migrated to the per-tool <tool>/variant_type layout: "
        + ", ".join(_PENDING) + ". Only fastVEP is available in this build.")


rule all:
    input:
        [str(p) for p in FASTVEP_TARGETS],
