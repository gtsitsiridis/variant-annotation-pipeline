"""AbSplice (tissue-specific splicing), gagneurlab/absplice.

AbSplice ships its own (heavy) snakemake workflow needing per-tissue SpliceMaps + SpliceAI
(+ GPU); we do NOT run it. Instead `config['absplice']['result']` points at a *precomputed*
AbSplice2 result — deeprvat's per-gene `*_max_preds.parquet` directory — and this rule
inner-joins the `AbSplice_DNA_max` / `AbSplice2_max` scores onto our variant set (SNV-only:
indels get no score). Only built when `absplice.enabled: true`.
"""


rule absplice:
    input:
        absplice_dir=config["absplice"]["result"],
        vep=OUT / "vep.parquet",
    output:
        parquet=OUT / "absplice" / "absplice.parquet",
    conda:
        "../../envs/parse.yaml"
    script:
        "../scripts/run_absplice.py"
