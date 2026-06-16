"""Tier 2 — AbSplice-DNA (tissue-specific splicing), gagneurlab/absplice.

AbSplice ships its own (heavy) snakemake workflow needing per-tissue SpliceMaps. We do NOT
vendor it: run AbSplice separately against your checkout, point `config['absplice']['result']`
at its output table, and this rule aggregates that to a per-variant AbSplice_DNA (max across
the requested tissues) keyed on the variant. Only built when `absplice.enabled: true`.
"""


rule absplice:
    input:
        result=config["absplice"]["result"],
    output:
        parquet=OUT / "absplice" / "absplice.parquet",
    params:
        tissues=",".join(config["absplice"].get("tissues", [])),
        per_tissue=config["absplice"].get("per_tissue", False),
    conda:
        "../../envs/parse.yaml"
    shell:
        "python workflow/scripts/run_absplice.py "
        "  --absplice-result {input.result} --tissues '{params.tissues}' "
        "  --per-tissue {params.per_tissue} --out {output.parquet}"
