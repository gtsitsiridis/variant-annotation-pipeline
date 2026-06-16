"""Tier 2 — AbSplice-DNA (tissue-specific splicing), gagneurlab/absplice.

AbSplice ships its own snakemake workflow + config (it needs SpliceMap reference data per
tissue). Rather than vendoring it, this rule shells out to an AbSplice run configured via
`config['absplice']['config']` and collects the per-variant, per-tissue AbSplice_DNA score
into a tidy parquet keyed on the variant. Only built when `absplice.enabled: true`.

Setup: clone gagneurlab/absplice, install its conda env, download SpliceMaps for the
requested tissues, and point `config['absplice']['config']` at its config.yaml.
"""


rule absplice:
    input:
        vcf=config["input_vcf"],
    output:
        parquet=OUT / "absplice" / "absplice.parquet",
    params:
        absplice_config=config["absplice"]["config"],
        tissues=",".join(config["absplice"]["tissues"]),
        outdir=lambda w, output: str(Path(output.parquet).parent),
    conda:
        "../../envs/absplice.yaml"
    shell:
        # Placeholder: drive the AbSplice workflow on {input.vcf}, then aggregate its
        # tissue outputs to a tidy parquet (chrom,start,ref,alt,tissue,AbSplice_DNA →
        # pivoted/maxed per variant). Wire to your AbSplice checkout.
        "python workflow/scripts/run_absplice.py "
        "  --vcf {input.vcf} --absplice-config {params.absplice_config} "
        "  --tissues {params.tissues} --out {output.parquet}"
