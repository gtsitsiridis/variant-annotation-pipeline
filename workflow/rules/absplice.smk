"""AbSplice2 — precomputed per-(variant, gene) splice scores (lookup, no compute).

Joins the precomputed AbSplice2 result (`additional.absplice.result`) onto fastVEP's filtered VCF
-> parts/absplice.parquet (variant_id, gene, AbSplice_DNA_max, AbSplice2_max). Gene-level, SNV-only;
combine joins it onto the fastVEP base on (variant_id, gene). Opt-in via `additional.absplice.enabled`.
"""

_ABSPLICE = config["additional"]["absplice"]


def _gb(x, default=16):
    """'16GB' / '16' -> 16 (int GiB); anything unparseable -> default."""
    try:
        return int(float(str(x).upper().replace("GB", "").strip()))
    except (TypeError, ValueError):
        return default


_ABSPLICE_MEM_GB = _gb(_ABSPLICE.get("memory_limit", "16GB"))
_ABSPLICE_THREADS = int(_ABSPLICE.get("threads", 4))


rule absplice:
    input:
        vcf=FASTVEP_FILTERED_VCF,
    output:
        parquet=PARTS / "absplice.parquet",
    params:
        result=_ABSPLICE["result"],       # precomputed AbSplice2 per-gene *_max_preds.parquet dir
        gtf=config["gtf"],                 # gencode: unversioned -> versioned gene_id map
        memory_limit=_ABSPLICE.get("memory_limit", "16GB"),
        threads=_ABSPLICE_THREADS,
    threads: _ABSPLICE_THREADS
    resources:
        mem_mb=_ABSPLICE_MEM_GB * 1024 + 8192,
        runtime=_ABSPLICE.get("runtime", 240),
    conda:
        "../../envs/parse.yaml"
    script:
        "../scripts/run_absplice.py"
