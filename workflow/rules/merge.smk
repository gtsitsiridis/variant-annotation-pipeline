"""Concatenate per-chromosome VEP parquets and merge in AbSplice → annotations.parquet."""


rule merge_vep:
    input:
        parquets=expand(OUT / "vep" / "{chrom}.vep.parquet", chrom=CHROMS),
    output:
        parquet=OUT / "vep.parquet",
    conda:
        "../../envs/parse.yaml"
    run:
        import polars as pl
        pl.concat([pl.scan_parquet(p) for p in input.parquets]).sink_parquet(output.parquet)


rule merge:
    """Final per-variant annotation table. AbSplice merged in only if enabled."""
    input:
        vep=OUT / "vep.parquet",
        absplice=(OUT / "absplice" / "absplice.parquet") if ABSPLICE else [],
    output:
        parquet=OUT / "annotations.parquet",
    params:
        with_absplice=ABSPLICE,
    conda:
        "../../envs/parse.yaml"
    script:
        "../scripts/merge_annotations.py"
