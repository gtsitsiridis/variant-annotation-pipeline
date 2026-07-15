"""Deep-tier selector — narrow the broad fastVEP set to the window the heavy tools score.

Reads the `small` track's basic annotation (Tier 0), keeps variants within `deep.window` bp of
a transcript (or inside one, or non-MODIFIER), and subsets the input VCF to those ids. The
result keeps the chrom_start_ref_alt ID convention, so the deep rules (vep/nmd/absplice) are
unchanged apart from reading this smaller VCF.
"""

DEEP_WINDOW = config["distance"]                 # the single pipeline cis window
ANALYSIS_SET = OUT / "deep" / "analysis_set.vcf.gz"


rule select_analysis_set_ids:
    """Pick variant_ids within `distance` of a gene from the fastVEP small partitions."""
    input:
        [FASTVEP_DIR / f"variant_type={vt}" for vt in ("SNV", "indel")],
    output:
        ids=OUT / "deep" / "analysis_set.ids.txt",
    params:
        fastvep_dir=str(FASTVEP_DIR),
        window=DEEP_WINDOW,
    conda:
        "../../envs/parse.yaml"
    script:
        "../scripts/select_analysis_set.py"


rule analysis_set:
    """Subset the input VCF to the selected ids -> analysis_set.vcf.gz (deep-tier input)."""
    input:
        vcf=config["input_vcf"],
        ids=OUT / "deep" / "analysis_set.ids.txt",
    output:
        vcf=ANALYSIS_SET,
        tbi=str(ANALYSIS_SET) + ".tbi",
    conda:
        "../../envs/vep.yaml"
    shell:
        "bcftools view -i 'ID=@{input.ids}' -Oz -o {output.vcf} {input.vcf} && "
        "tabix -p vcf {output.vcf}"
