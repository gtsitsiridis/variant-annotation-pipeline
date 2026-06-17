"""NMD-Scanner (gagneurlab/NMD-Scanner) — NMD escape prediction for PTC variants.

Per transcript×variant: detects premature termination codons and evaluates the five
canonical NMD escape rules (last-exon, 50nt-penultimate, long-exon, start-proximal,
single-exon) → `nmd_escape`. Refines the pLoF signal beyond LOFTEE (an NMD-escaping
stop_gained keeps a truncated protein; an NMD-triggering one degrades the transcript).

Scattered per chunk, then normalized to the shared variant key (+ transcript) for the
merge. We request CSV output: NMD-Scanner's parquet writer chokes on a mixed-type column
(`ref_all_stop_codons`), and we only keep a few boolean columns downstream anyway.
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
        csv=OUT / "nmd" / "chunk_{chunk}.nmd.csv",
    params:
        reassign=_REASSIGN,
    conda:
        "../../envs/nmd.yaml"
    shell:
        "nmd-scanner --vcf {input.vcf} --gtf {input.gtf} --fasta {input.fasta} "
        "  --output {output.csv} {params.reassign}"


def _nmd_chunk_csvs(wildcards):
    return expand(str(OUT / "nmd" / "chunk_{chunk}.nmd.csv"), chunk=chunk_ids())


rule merge_nmd:
    """Concat per-chunk NMD CSVs, normalize keys to (chrom,start,ref,alt,Feature)."""
    input:
        csvs=_nmd_chunk_csvs,
    output:
        parquet=OUT / "nmd.parquet",
    run:
        import polars as pl
        frames = []
        for p in input.csvs:
            # Project only the columns we keep — pushdown skips NMD-Scanner's messy
            # mixed-type columns (e.g. ref_all_stop_codons), so the wide CSV scans cleanly.
            lf = pl.scan_csv(p)
            have = set(lf.collect_schema().names())
            keep = [c for c in NMD_COLS if c in have]
            frames.append(lf.select(
                pl.col("chromosome").alias("chrom"),
                # NMD-Scanner's start_variant is already 0-based (POS-1) — the same coordinate
                # as our VEP variant key — so use it directly. Verified: VCF POS 69869 ->
                # ID chr1_69868_T_A -> VEP start 69868 == NMD start_variant 69868.
                pl.col("start_variant").cast(pl.Int64, strict=False).alias("start"),
                pl.col("ref"), pl.col("alt"),
                pl.col("transcript_id").alias("Feature"),
                # Cast to a stable Boolean so empty chunks don't break the concat schema.
                *[(pl.col(c).cast(pl.String).str.to_lowercase() == "true").alias(c) for c in keep],
            ))
        pl.concat(frames).sink_parquet(output.parquet)
