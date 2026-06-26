"""Tier 3 — ENCODE-rE2G enhancer -> gene annotation (hg38, precomputed).

Per variant: which ENCODE-rE2G enhancer it falls in, the predicted target gene, the rE2G
score, and the variant->TSS distance. STANDALONE parquet (grain = variant x target_gene),
joined by the consumer on variant_id (NOT folded into the transcript-level annotations.parquet
— different grain). Opt-in via config `e2g.enabled`. Reads the on-disk E2G export
(predictions dir + enhancers.parquet); no liftOver (already hg38).
"""

_E2G = config.get("e2g", {})


rule e2g:
    input:
        vcf=config["input_vcf"],
    output:
        parquet=OUT / "e2g" / "e2g.parquet",
    params:
        predictions=_E2G["predictions"],          # dir of per-chrom prediction parquets
        enhancers=_E2G["enhancers"],               # enhancers.parquet (id -> chrom/start/end/class)
        gtf=config["gtf"],                         # gencode: unversioned -> versioned gene_id map
        model=_E2G.get("model", "ENCODE-rE2G"),    # vs scE2G
        score_threshold=_E2G.get("score_threshold", 0.5),
        classes=_E2G.get("classes", ["genic", "intergenic"]),   # distal; drop promoter
        cell_types=_E2G.get("cell_types", None),   # None = all biosamples; or a list to restrict
        memory_limit=_E2G.get("memory_limit", "32GB"),
        threads=_E2G.get("threads", 4),
    conda:
        "../../envs/e2g.yaml"
    script:
        "../scripts/annotate_enhancers.py"
