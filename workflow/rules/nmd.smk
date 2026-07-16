"""NMD-Scanner (gagneurlab/NMD-Scanner) — NMD escape prediction for PTC variants.

Per transcript×variant: detects premature termination codons and evaluates the five canonical NMD
escape rules (last-exon, 50nt-penultimate, long-exon, start-proximal, single-exon) -> `nmd_escape`.
Refines the pLoF signal beyond LOFTEE (an NMD-escaping stop_gained keeps a truncated protein; an
NMD-triggering one degrades the transcript).

Runs on fastVEP's filtered VCF -> parts/nmd.parquet (variant_id, transcript + NMD flags); combine
joins it onto the base on (variant_id, transcript). Opt-in via `additional.nmd.enabled`. We request
CSV output: NMD-Scanner's parquet writer chokes on a mixed-type column (`ref_all_stop_codons`), and
we only keep a few boolean columns downstream. (Single job on the filtered subset; chunk like VEP
if it becomes a straggler.)
"""

_NMD = config["additional"]["nmd"]
_REASSIGN = "--reassign_exons" if _NMD.get("reassign_exons", False) else ""

# Boolean columns kept from NMD-Scanner's wide output (kept if present).
NMD_COLS = [
    "nmd_escape", "alt_is_premature", "start_loss", "stop_loss",
    "nmd_last_exon_rule", "nmd_50nt_penultimate_rule", "nmd_long_exon_rule",
    "nmd_start_proximal_rule", "nmd_single_exon_rule",
]


rule nmd:
    """NMD-Scanner on the filtered VCF -> wide CSV (temp)."""
    input:
        vcf=FASTVEP_FILTERED_VCF,
        gtf=config["gtf"],
        fasta=config["fasta"],
    output:
        csv=temp(OUT / "nmd" / "nmd.csv"),
    params:
        reassign=_REASSIGN,
    resources:
        mem_mb=_NMD.get("mem_mb", 32000),
        runtime=_NMD.get("runtime", 720),
    conda:
        "../../envs/nmd.yaml"
    shell:
        "nmd-scanner --vcf {input.vcf} --gtf {input.gtf} --fasta {input.fasta} "
        "  --output {output.csv} {params.reassign}"


rule parse_nmd:
    """NMD CSV -> parts/nmd.parquet, keyed on (variant_id, transcript) for the combine join."""
    input:
        csv=OUT / "nmd" / "nmd.csv",
    output:
        parquet=PARTS / "nmd.parquet",
    run:
        import polars as pl
        lf = pl.scan_csv(input.csv)             # project only the columns we keep (skip messy wide cols)
        have = set(lf.collect_schema().names())
        keep = [c for c in NMD_COLS if c in have]
        (lf.select(
            # variant_id = chrom_start_ref_alt (start_variant is already 0-based POS-1, == our key).
            variant_id=pl.concat_str([
                pl.col("chromosome"),
                pl.col("start_variant").cast(pl.Int64, strict=False).cast(pl.String),
                pl.col("ref"), pl.col("alt")], separator="_"),
            transcript=pl.col("transcript_id"),
            *[(pl.col(c).cast(pl.String).str.to_lowercase() == "true").alias(c) for c in keep],
        ).sink_parquet(output.parquet))
