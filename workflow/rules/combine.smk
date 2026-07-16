"""Combine step — the fastVEP base + enabled additional tools -> annotations.parquet.

LEFT-joins the per-method parts (fastvep base + whichever of vep/e2g/absplice are enabled) into ONE
table and writes it hive-partitioned by variant_type / gene_type / canonical. Always runs (with no
additional tools it is just the filtered + partitioned fastVEP base). Reads the module-level
VEP/E2G/ABSPLICE flags + FASTVEP_* from fastvep.smk / the Snakefile.
"""

ANNOTATIONS_DIR = OUT / "annotations.parquet"
# variant_type partitions that exist given the configured tracks (static; gene_type/canonical nest below).
COMBINE_VTYPES = ["SNV", "indel"] + (["SV"] if "sv" in TRACKS else [])


def _combine_inputs(wc):
    d = {"fastvep": str(FASTVEP_PART)}
    if VEP:
        d["vep"] = str(PARTS / "vep.parquet")
    if E2G:
        d["e2g"] = str(PARTS / "e2g.parquet")
    if ABSPLICE:
        d["absplice"] = str(PARTS / "absplice.parquet")
    return d


rule combine:
    input:
        unpack(_combine_inputs),
    output:
        [directory(ANNOTATIONS_DIR / f"variant_type={vt}") for vt in COMBINE_VTYPES],
    params:
        out_dir=str(ANNOTATIONS_DIR),
        protein_coding_only=FASTVEP_PROTEIN_CODING_ONLY,
        canonical_only=FASTVEP_CANONICAL_ONLY,
        memory_limit=_FASTVEP.get("parse_memory_limit", "64GB"),
        threads=_FASTVEP.get("parse_threads", 6),
    conda:
        "../../envs/parse.yaml"
    script:
        "../scripts/combine_annotations.py"
