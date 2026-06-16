"""NMD-Scanner (gagneurlab/NMD-Scanner) — NMD escape prediction for PTC variants.

Per transcript×variant: detects premature termination codons and evaluates the five
canonical NMD escape rules (last-exon, 50nt-penultimate, long-exon, start-proximal,
single-exon) → `nmd_escape`. Refines the pLoF signal beyond LOFTEE (an NMD-escaping
stop_gained keeps a truncated protein; an NMD-triggering one degrades the transcript).

Scattered per chromosome, then normalized to the shared variant key (+ transcript) for the
merge. NMD-Scanner writes parquet natively (needs pyarrow, in envs/nmd.yaml).
"""

_REASSIGN = "--reassign_exons" if config.get("nmd", {}).get("reassign_exons", False) else ""

# Columns we keep from NMD-Scanner's wide output (kept if present).
NMD_COLS = [
    "nmd_escape", "alt_is_premature", "start_loss", "stop_loss",
    "nmd_last_exon_rule", "nmd_50nt_penultimate_rule", "nmd_long_exon_rule",
    "nmd_start_proximal_rule", "nmd_single_exon_rule",
]


rule nmd:
    input:
        vcf=CHUNKS / "chunk_{chunk}.vcf.gz",
        gtf=config["gtf"],
        fasta=config["fasta"],
    output:
        parquet=OUT / "nmd" / "chunk_{chunk}.nmd.parquet",
    params:
        reassign=_REASSIGN,
    conda:
        "../../envs/nmd.yaml"
    shell:
        "nmd-scanner --vcf {input.vcf} --gtf {input.gtf} --fasta {input.fasta} "
        "  --output {output.parquet} {params.reassign}"


def _nmd_chunk_parquets(wildcards):
    return expand(str(OUT / "nmd" / "chunk_{chunk}.nmd.parquet"), chunk=chunk_ids())


rule merge_nmd:
    """Concat per-chunk NMD output, normalize keys to (chrom,start,ref,alt,Feature)."""
    input:
        parquets=_nmd_chunk_parquets,
    output:
        parquet=OUT / "nmd.parquet",
    run:
        import polars as pl
        frames = []
        for p in input.parquets:
            lf = pl.scan_parquet(p)
            have = set(lf.collect_schema().names())
            keep = [c for c in NMD_COLS if c in have]
            frames.append(lf.select(
                pl.col("chromosome").alias("chrom"),
                # NMD start_variant is the VCF POS (1-based); our variant key `start` is
                # 0-based (pos = start+1). VERIFY indel representation when smoke-testing.
                (pl.col("start_variant") - 1).alias("start"),
                pl.col("ref"), pl.col("alt"),
                pl.col("transcript_id").alias("Feature"),
                *[pl.col(c) for c in keep],
            ))
        pl.concat(frames).sink_parquet(output.parquet)
