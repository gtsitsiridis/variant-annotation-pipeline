"""fastVEP funnel — the always-on driver.

fastVEP annotates the WHOLE variant set (small + optional SV) cheaply (only gencode GFF3 + FASTA),
then serves two roles:
  1. the base annotation — the CSQ is exploded + joined to transcript_metadata into ONE flat table
     `parts/fastvep.parquet` (variant × transcript, all variant_types), consumed by combine.
  2. the funnel — the small variants are filtered (distance + protein_coding_only + canonical_only)
     to `fastvep/small.vcf.gz`, the VCF that every `additional` tool (vep/e2g/absplice) re-annotates.

Tracks: `small` (config `input_vcf`, ID=chrom_start_ref_alt) and, if set, `sv` (config `sv_vcf`,
ID=sv_id). SVs are fastVEP-only — they enter parts/fastvep.parquet but never the filtered VCF.
"""

_FASTVEP = config["fastvep"]
FASTVEP_DISTANCE = _FASTVEP["distance"]                      # cis window: fastVEP --distance + filter
FASTVEP_PROTEIN_CODING_ONLY = _FASTVEP.get("protein_coding_only", False)
FASTVEP_CANONICAL_ONLY = _FASTVEP.get("canonical_only", False)
FASTVEP_INCLUDE_SV = _FASTVEP.get("include_sv", True)

# track -> input VCF. `small` is mandatory; `sv` only if sv_vcf is set AND include_sv.
TRACK_VCF = {"small": config["input_vcf"]}
if config.get("sv_vcf") and FASTVEP_INCLUDE_SV:
    TRACK_VCF["sv"] = config["sv_vcf"]
TRACKS = list(TRACK_VCF)

# Huge temporary CSQ VCFs (bgzipped, temp()) — keep them off the network output dir.
FASTVEP_SCRATCH = _FASTVEP.get("scratch", str(OUT / "fastvep"))

PARTS = OUT / "parts"
FASTVEP_PART = PARTS / "fastvep.parquet"                    # flat transcript table (base for filter + combine)
FASTVEP_IDS = OUT / "fastvep" / "small.ids.txt"
FASTVEP_FILTERED_VCF = OUT / "fastvep" / "small.vcf.gz"     # THE fastVEP output VCF for additional tools


wildcard_constraints:
    track="|".join(TRACKS),


rule transcript_metadata:
    """gencode GTF -> transcript metadata (gene_type/symbol/biotype/canonical/tsl) for the CSQ join."""
    input:
        gtf=config["gtf"],
    output:
        parquet=OUT / "transcript_metadata.parquet",
    conda:
        "../../envs/parse.yaml"
    script:
        "../scripts/build_transcript_metadata.py"


rule fastvep:
    """fastVEP annotate -> bgzipped CSQ VCF (all transcripts per variant within --distance), temp."""
    input:
        vcf=lambda wc: TRACK_VCF[wc.track],
        gff3=config["gff3"],
        fasta=config["fasta"],
    output:
        vcf=temp(f"{FASTVEP_SCRATCH}/{{track}}.csq.vcf.gz"),
    params:
        distance=FASTVEP_DISTANCE,
    threads: 4
    conda:
        "../../envs/fastvep.yaml"
    shell:
        "fastvep annotate -i {input.vcf} --gff3 {input.gff3} --fasta {input.fasta} "
        "  --output-format vcf --everything --hgvs --canonical --symbol "
        "  --distance {params.distance} -o /dev/stdout | bgzip -@ {threads} > {output.vcf}"


def _parse_inputs(wc):
    """small CSQ (always) + metadata (+ sv CSQ when the SV track is configured)."""
    d = {"small": f"{FASTVEP_SCRATCH}/small.csq.vcf.gz",
         "metadata": str(OUT / "transcript_metadata.parquet")}
    if "sv" in TRACKS:
        d["sv"] = f"{FASTVEP_SCRATCH}/sv.csq.vcf.gz"
    return d


rule parse_fastvep:
    """Explode the CSQ VCF(s) -> ONE flat parts/fastvep.parquet (variant × transcript, all variant_types)."""
    input:
        unpack(_parse_inputs),
    output:
        parquet=FASTVEP_PART,
    params:
        memory_limit=_FASTVEP.get("parse_memory_limit", "32GB"),
        threads=_FASTVEP.get("parse_threads", 4),
    conda:
        "../../envs/parse.yaml"
    script:
        "../scripts/parse_fastvep.py"


rule fastvep_filter_ids:
    """Filter the small transcript rows (distance + protein_coding_only + canonical_only) -> variant ids."""
    input:
        fastvep=FASTVEP_PART,
    output:
        ids=FASTVEP_IDS,
    params:
        window=FASTVEP_DISTANCE,
        protein_coding_only=FASTVEP_PROTEIN_CODING_ONLY,
        canonical_only=FASTVEP_CANONICAL_ONLY,
    conda:
        "../../envs/parse.yaml"
    script:
        "../scripts/select_analysis_set.py"


rule fastvep_filtered_vcf:
    """Subset input_vcf to the passing ids -> fastvep/small.vcf.gz (input to all additional tools)."""
    input:
        vcf=config["input_vcf"],
        ids=FASTVEP_IDS,
    output:
        vcf=FASTVEP_FILTERED_VCF,
        tbi=str(FASTVEP_FILTERED_VCF) + ".tbi",
    conda:
        "../../envs/vep.yaml"
    shell:
        "bcftools view -i 'ID=@{input.ids}' -Oz -o {output.vcf} {input.vcf} && "
        "tabix -p vcf {output.vcf}"
