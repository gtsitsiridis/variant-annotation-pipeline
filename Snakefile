"""variant_annotation — dataset-agnostic functional annotation of a variant set.

fastVEP is the always-on DRIVER: it annotates the whole set (small variants ID=chrom_start_ref_alt;
optional SVs ID=sv_id), the CSQ is exploded + joined to transcript_metadata into ONE flat base table
(parts/fastvep.parquet), and the small variants are filtered (distance + protein_coding_only +
canonical_only) to `fastvep/small.vcf.gz`. Every enabled `additional` tool re-annotates that filtered
VCF, and combine LEFT-joins them all into ONE `annotations.parquet`, hive-partitioned by
variant_type / gene_type / canonical.

  fastVEP   ALWAYS on (gff3 + fasta). Base annotation + the funnel; the only SV-capable tool.
  vep       additional (additional.vep.enabled). LOFTEE/CADD/SpliceAI/PrimateAI/AlphaMissense — transcript-level.
  nmd       additional (additional.nmd.enabled). NMD-Scanner escape prediction for PTCs — transcript-level.
  e2g       additional (additional.e2g.enabled). ENCODE-rE2G enhancer -> target gene — gene-level.
  absplice  additional (additional.absplice.enabled). Precomputed AbSplice2 lookup — gene-level, SNV-only.

Configure via config.yaml (or --configfile), or pass config as a Snakemake `module`.
"""
from pathlib import Path

configfile: "config.yaml"

OUT = Path(config["output_dir"])

# ── additional-tool enable flags (each re-annotates fastVEP's filtered VCF) ───────────
_ADD = config.get("additional", {})
VEP = _ADD.get("vep", {}).get("enabled", False)
NMD = _ADD.get("nmd", {}).get("enabled", False)
E2G = _ADD.get("e2g", {}).get("enabled", False)
ABSPLICE = _ADD.get("absplice", {}).get("enabled", False)

# SV annotation is fastVEP-only. Reject include_sv on any additional tool so the flag never lies.
for _t in ("vep", "nmd", "e2g", "absplice"):
    if _ADD.get(_t, {}).get("include_sv", False):
        raise ValueError(f"additional.{_t}.include_sv is not supported — SV annotation is fastVEP-only.")

wildcard_constraints:
    chunk=r"\d+",

include: "workflow/rules/fastvep.smk"        # base + funnel: parts/fastvep.parquet + fastvep/small.vcf.gz
if VEP:
    include: "workflow/rules/vep.smk"        # -> parts/vep.parquet
if NMD:
    include: "workflow/rules/nmd.smk"        # -> parts/nmd.parquet
if E2G:
    include: "workflow/rules/e2g.smk"        # -> parts/e2g.parquet
if ABSPLICE:
    include: "workflow/rules/absplice.smk"   # -> parts/absplice.parquet
include: "workflow/rules/combine.smk"        # -> annotations.parquet (always; base + enabled tools)


rule all:
    input:
        [str(ANNOTATIONS_DIR / f"variant_type={vt}") for vt in COMBINE_VTYPES],
