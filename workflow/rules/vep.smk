"""Tier 1 — Ensembl VEP + plugins, scattered per chromosome.

NOTE on gene namespace: to match a specific gencode release (gtex-benchmark = v34) either
use a VEP cache for the matching Ensembl release, or force the custom GTF with
`--gtf {gtf} --fasta {fasta}` (uncomment below). The cache route is shown by default
because the plugins (esp. LOFTEE) are best-supported against the cache.
"""

VEP_PD = config["plugin_data"]


def _plugin_args() -> str:
    """Assemble the --plugin flags from the configured reference data."""
    pd = VEP_PD
    plugins = [
        f"--plugin CADD,{pd['cadd_snv']},{pd['cadd_indel']}",
        f"--plugin SpliceAI,snv={pd['spliceai_snv']},indel={pd['spliceai_indel']}",
        f"--plugin PrimateAI,{pd['primateai']}",
        f"--plugin AlphaMissense,file={pd['alphamissense']}",
        f"--plugin Condel,{pd['condel_config']},s,2",
        # LOFTEE (loftee dir must be on --dir_plugins or PERL5LIB)
        (f"--plugin LoF,loftee_path:{pd['loftee_path']},"
         f"human_ancestor_fa:{pd['loftee_human_ancestor']},"
         f"conservation_file:{pd['loftee_conservation']},"
         f"gerp_bigwig:{pd['loftee_gerp']}"),
    ]
    return " ".join(plugins)


rule split_vcf:
    """One bgzipped VCF per chromosome (keeps each VEP job tractable)."""
    input:
        vcf=config["input_vcf"],
    output:
        vcf=OUT / "vcf" / "{chrom}.vcf.gz",
    conda:
        "../../envs/vep.yaml"
    shell:
        "bcftools view -r {wildcards.chrom} -Oz -o {output.vcf} {input.vcf} && "
        "tabix -p vcf {output.vcf}"


rule vep:
    """VEP + plugins → tab output. Uploaded_variation = chrom_start_ref_alt (join key)."""
    input:
        vcf=OUT / "vcf" / "{chrom}.vcf.gz",
        fasta=config["fasta"],
    output:
        tab=OUT / "vep" / "{chrom}.vep.tsv.gz",
    params:
        cache_dir=config["vep"]["cache_dir"],
        plugin_dir=config["vep"]["plugin_dir"],
        cache_version=config["vep"]["cache_version"],
        assembly=config["assembly"],
        fork=config["vep"]["fork"],
        plugins=_plugin_args,
        loftee_path=VEP_PD["loftee_path"],
    conda:
        "../../envs/vep.yaml"
    shell:
        # PERL5LIB must include LOFTEE for its modules to load.
        "PERL5LIB={params.loftee_path}:$PERL5LIB "
        "vep --offline --cache --dir_cache {params.cache_dir} "
        "    --dir_plugins {params.plugin_dir} --cache_version {params.cache_version} "
        "    --fasta {input.fasta} --assembly {params.assembly} --fork {params.fork} "
        "    --everything --force_overwrite --tab --compress_output bgzip "
        "    --input_file {input.vcf} --output_file {output.tab} "
        "    {params.plugins}"
        # To force a gencode build instead of the cache, replace --cache/--dir_cache with:
        #   --gtf {config[gtf]}   (requires bgzip+tabix'd, sorted GTF)


rule parse_vep:
    """VEP tab → tidy parquet (one row per variant×transcript), keyed on the variant."""
    input:
        tab=OUT / "vep" / "{chrom}.vep.tsv.gz",
    output:
        parquet=OUT / "vep" / "{chrom}.vep.parquet",
    conda:
        "../../envs/parse.yaml"
    script:
        "../scripts/parse_vep.py"
