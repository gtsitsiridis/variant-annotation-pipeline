"""ENCODE-rE2G enhancer -> gene annotation (hg38, precomputed).

Per (variant, target gene): which ENCODE-rE2G enhancer the variant falls in, the predicted target
gene, the rE2G score, and the variant->TSS distance. Gene-level — combine joins it onto the fastVEP
base on (variant_id, gene). Reads fastVEP's filtered VCF (fastvep/small.vcf.gz) like the other
additional tools; opt-in via `additional.e2g.enabled`. No liftOver (the export is already hg38).
"""

_E2G = config["additional"]["e2g"]


def _gb(x, default=32):
    """'32GB' / '32' -> 32 (int GiB); anything unparseable -> default."""
    try:
        return int(float(str(x).upper().replace("GB", "").strip()))
    except (TypeError, ValueError):
        return default


_E2G_MEM_GB = _gb(_E2G.get("memory_limit", "32GB"))
_E2G_THREADS = int(_E2G.get("threads", 4))


rule e2g:
    input:
        vcf=FASTVEP_FILTERED_VCF,
    output:
        parquet=PARTS / "e2g.parquet",
    params:
        predictions=_E2G["predictions"],          # dir of per-chrom prediction parquets
        enhancers=_E2G["enhancers"],               # enhancers.parquet (id -> chrom/start/end/class)
        gtf=config["gtf"],                         # gencode: unversioned -> versioned gene_id map
        model=_E2G.get("model", "ENCODE-rE2G"),    # vs scE2G
        score_threshold=_E2G.get("score_threshold", 0.5),
        classes=_E2G.get("classes", ["genic", "intergenic"]),   # distal; drop promoter
        cell_types=_E2G.get("cell_types", None),   # None = all biosamples; or a list to restrict
        distance_to_tss=_E2G.get("distance_to_tss", None),      # None = no cap on variant->TSS distance
        memory_limit=_E2G.get("memory_limit", "32GB"),
        threads=_E2G_THREADS,
    threads: _E2G_THREADS
    resources:
        # SLURM mem must exceed DuckDB's memory_limit (it only spills to duckdb_tmp once it HITS the
        # limit) + headroom for python/arrow; else the job OOM-kills before spilling.
        mem_mb=_E2G_MEM_GB * 1024 + 8192,
        runtime=_E2G.get("runtime", 240),
    conda:
        "../../envs/e2g.yaml"
    script:
        "../scripts/annotate_enhancers.py"
