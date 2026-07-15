"""variant_annotation — dataset-agnostic functional annotation of a variant set.

Per-tool annotation over a generic VCF (small variants ID=chrom_start_ref_alt; SVs ID=sv_id),
under a distance-keyed layout: OUT/distance_<distance>/<tool>.parquet/variant_type=<vtype>/,
each hive-partitioned by variant_type and keyed on `variant_id`.

  fastVEP   ALWAYS on (only gff3+fasta). Broad per-transcript consequence/distance over the
            whole set; the only SV-capable tool (fastvep.include_sv). Base for the `distance`
            selection used by the other tools.  -> distance_<d>/fastvep.parquet/
  VEP       opt-in (vep.enabled). LOFTEE/CADD/SpliceAI/PrimateAI/AlphaMissense on the small set
            selected within `distance` of a gene (chunked).  -> distance_<d>/vep.parquet/
  NMD / AbSplice / E2G   opt-in — PENDING migration to this layout (enabling one errors).

A single top-level `distance` drives fastVEP `--distance`, the selection window, and VEP
`--distance`. Configure via config.yaml (or --configfile), or pass config as a Snakemake `module`.
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

# fastVEP + VEP are migrated to the per-tool distance_<d>/<tool>.parquet/variant_type layout.
# VEP selects its (small) analysis set from the fastVEP `distance` annotation, chunks + runs VEP,
# and writes distance_<d>/vep.parquet/variant_type={SNV,indel}/ keyed on variant_id.
if VEP:
    include: "workflow/rules/select.smk"     # ANALYSIS_SET from the fastVEP distance partitions
    include: "workflow/rules/vep.smk"        # defines VEP_DIR + the vep_parquet target
if E2G:
    include: "workflow/rules/e2g.smk"        # standalone (reads input_vcf + ENCODE-rE2G tables); defines E2G_DIR

# nmd / absplice are still pending migration to the new layout; enabling one fails loudly
# rather than running the retired deep-tier machinery.
_PENDING = [t for t, on in (("nmd", NMD), ("absplice", ABSPLICE)) if on]
if _PENDING:
    raise ValueError(
        "tools not yet migrated to the per-tool <tool>/variant_type layout: "
        + ", ".join(_PENDING) + ". Available: fastVEP, VEP, e2g.")


def _targets():
    t = list(FASTVEP_TARGETS)
    if VEP:
        t += [VEP_DIR / f"variant_type={vt}" for vt in ("SNV", "indel")]
    if E2G:
        t += [E2G_DIR / f"variant_type={vt}" for vt in ("SNV", "indel")]
    return [str(x) for x in t]


rule all:
    input:
        _targets(),
