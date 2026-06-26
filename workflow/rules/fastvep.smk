"""Tier 0 — fastVEP: fast per-transcript consequence/distance over the WHOLE variant set.

Cheap (only the gencode GFF3 + FASTA, no plugin reference data) so it always runs: it is both
the broad annotation layer and the basis for selecting the deep-tier window. Runs per input
track — `small` (config `input_vcf`, ID=chrom_start_ref_alt) and, if set, `sv` (config
`sv_vcf`, ID=sv_id). fastVEP is ID-agnostic (it round-trips the VCF ID as the join key), so
one rule annotates both; only the `small` track feeds the deep tiers downstream.

Per-track outputs land under OUT/<track>/basic_annotations.parquet so the {track} wildcard
stays clean.
"""

# track -> input VCF. `small` is mandatory; `sv` only if sv_vcf is configured.
TRACK_VCF = {"small": config["input_vcf"]}
if config.get("sv_vcf"):
    TRACK_VCF["sv"] = config["sv_vcf"]
TRACKS = list(TRACK_VCF)

_FASTVEP = config.get("fastvep", {})
FASTVEP_DISTANCE = _FASTVEP.get("distance", 5000)
# Huge temporary CSQ VCFs (the small track is ~20 GB) — keep them off the network output dir.
FASTVEP_SCRATCH = _FASTVEP.get("scratch", str(OUT / "fastvep"))

# The `small` track's basic annotation is the input to the deep-tier selector.
BASIC_SMALL = OUT / "small" / "basic_annotations.parquet"


wildcard_constraints:
    track="|".join(TRACKS),


rule transcript_metadata:
    """gencode GTF -> transcript metadata (symbol/biotype/canonical/tsl) for the CSQ join."""
    input:
        gtf=config["gtf"],
    output:
        parquet=OUT / "transcript_metadata.parquet",
    conda:
        "../../envs/parse.yaml"
    script:
        "../scripts/build_transcript_metadata.py"


rule fastvep:
    """fastVEP annotate -> CSQ VCF (all transcripts per variant within --distance)."""
    input:
        vcf=lambda wc: TRACK_VCF[wc.track],
        gff3=config["gff3"],
        fasta=config["fasta"],
    output:
        vcf=temp(f"{FASTVEP_SCRATCH}/{{track}}.csq.vcf"),
    params:
        distance=FASTVEP_DISTANCE,
    conda:
        "../../envs/fastvep.yaml"
    shell:
        "fastvep annotate -i {input.vcf} --gff3 {input.gff3} --fasta {input.fasta} "
        "  --output-format vcf --everything --hgvs --canonical --symbol "
        "  --distance {params.distance} -o {output.vcf}"


rule parse_fastvep:
    """Explode CSQ -> tidy transcript-level parquet keyed on variant_id (+ metadata join)."""
    input:
        vcf=f"{FASTVEP_SCRATCH}/{{track}}.csq.vcf",
        metadata=OUT / "transcript_metadata.parquet",
    output:
        parquet=OUT / "{track}" / "basic_annotations.parquet",
    params:
        memory_limit=_FASTVEP.get("parse_memory_limit", "32GB"),
        threads=_FASTVEP.get("parse_threads", 4),
    conda:
        "../../envs/parse.yaml"
    script:
        "../scripts/parse_fastvep.py"
