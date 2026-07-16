"""Tier 3 — ENCODE-rE2G enhancer -> gene annotation (hg38, precomputed).

Per variant: which ENCODE-rE2G enhancer it falls in, the predicted target gene, the rE2G
score, and the variant->TSS distance. STANDALONE parquet (grain = variant x target_gene),
joined by the consumer on variant_id (NOT folded into the transcript-level annotations.parquet
— different grain). Opt-in via config `e2g.enabled`. Reads the on-disk E2G export
(predictions dir + enhancers.parquet); no liftOver (already hg38).
"""

_E2G = config.get("e2g", {})
E2G_DIR = OUT / f"distance_{config['distance']}" / "e2g.parquet"   # hive by variant_type (small-only)


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
        vcf=config["input_vcf"],
    output:
        [directory(E2G_DIR / f"variant_type={vt}") for vt in ("SNV", "indel")],
    params:
        out_dir=str(E2G_DIR),
        predictions=_E2G["predictions"],          # dir of per-chrom prediction parquets
        enhancers=_E2G["enhancers"],               # enhancers.parquet (id -> chrom/start/end/class)
        gtf=config["gtf"],                         # gencode: unversioned -> versioned gene_id map
        model=_E2G.get("model", "ENCODE-rE2G"),    # vs scE2G
        score_threshold=_E2G.get("score_threshold", 0.5),
        classes=_E2G.get("classes", ["genic", "intergenic"]),   # distal; drop promoter
        cell_types=_E2G.get("cell_types", None),   # None = all biosamples; or a list to restrict
        memory_limit=_E2G.get("memory_limit", "32GB"),
        threads=_E2G_THREADS,
    threads: _E2G_THREADS
    resources:
        # SLURM mem must exceed DuckDB's memory_limit (it only spills to duckdb_tmp once it
        # HITS the limit) + headroom for python/arrow; else the job OOM-kills before spilling.
        mem_mb=_E2G_MEM_GB * 1024 + 8192,
        runtime=_E2G.get("runtime", 240),
    conda:
        "../../envs/e2g.yaml"
    script:
        "../scripts/annotate_enhancers.py"
