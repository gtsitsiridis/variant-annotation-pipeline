"""Tier 0 — fastVEP: fast per-transcript consequence/distance over the WHOLE variant set.

Cheap (only the gencode GFF3 + FASTA, no plugin reference data) so it always runs: it is both
the broad annotation layer and the basis for selecting the deep-tier window. Runs per input
track — `small` (config `input_vcf`, ID=chrom_start_ref_alt) and, if set, `sv` (config
`sv_vcf`, ID=sv_id). fastVEP is ID-agnostic (it round-trips the VCF ID as the join key), so
one rule annotates both; only the `small` track feeds the deep tiers downstream.

Per-track outputs land under OUT/<track>/basic_annotations.parquet so the {track} wildcard
stays clean.
"""

# track -> input VCF. `small` is mandatory; `sv` only if sv_vcf is set AND fastvep.include_sv.
_FASTVEP = config.get("fastvep", {})
FASTVEP_DISTANCE = config["distance"]            # the single pipeline cis window (also prefixes the output dir)
FASTVEP_INCLUDE_SV = _FASTVEP.get("include_sv", True)
TRACK_VCF = {"small": config["input_vcf"]}
if config.get("sv_vcf") and FASTVEP_INCLUDE_SV:
    TRACK_VCF["sv"] = config["sv_vcf"]
TRACKS = list(TRACK_VCF)

# Huge temporary CSQ VCFs (the small track is ~20 GB) — keep them off the network output dir.
FASTVEP_SCRATCH = _FASTVEP.get("scratch", str(OUT / "fastvep"))

# Output: OUT/distance_<distance>/fastvep.parquet/variant_type={SNV,indel,SV}/<track>.parquet
# (hive-partitioned). `small` -> {SNV,indel} (split by ref/alt in parse); `sv` -> {SV}. This is the
# base annotation table every downstream tool + consumer reads
# (distance_<d>/fastvep.parquet/**/*.parquet, hive_partitioning=true).
FASTVEP_DIR = OUT / f"distance_{FASTVEP_DISTANCE}" / "fastvep.parquet"
FASTVEP_PARTS = {"small": ["SNV", "indel"], "sv": ["SV"]}
# variant_type partition dirs that exist given the configured tracks (Tier-0 targets).
FASTVEP_TARGETS = [FASTVEP_DIR / f"variant_type={vt}" for tr in TRACKS for vt in FASTVEP_PARTS[tr]]


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


rule parse_fastvep_small:
    """Explode the small-track CSQ -> fastvep/variant_type={SNV,indel}/ (split by ref/alt)."""
    input:
        vcf=f"{FASTVEP_SCRATCH}/small.csq.vcf",
        metadata=OUT / "transcript_metadata.parquet",
    output:
        directory(FASTVEP_DIR / "variant_type=SNV"),
        directory(FASTVEP_DIR / "variant_type=indel"),
    params:
        out_dir=str(FASTVEP_DIR),
        track="small",
        memory_limit=_FASTVEP.get("parse_memory_limit", "32GB"),
        threads=_FASTVEP.get("parse_threads", 4),
    conda:
        "../../envs/parse.yaml"
    script:
        "../scripts/parse_fastvep.py"


if "sv" in TRACKS:

    rule parse_fastvep_sv:
        """Explode the SV-track CSQ -> fastvep/variant_type=SV/."""
        input:
            vcf=f"{FASTVEP_SCRATCH}/sv.csq.vcf",
            metadata=OUT / "transcript_metadata.parquet",
        output:
            directory(FASTVEP_DIR / "variant_type=SV"),
        params:
            out_dir=str(FASTVEP_DIR),
            track="sv",
            memory_limit=_FASTVEP.get("parse_memory_limit", "32GB"),
            threads=_FASTVEP.get("parse_threads", 4),
        conda:
            "../../envs/parse.yaml"
        script:
            "../scripts/parse_fastvep.py"
