"""Concatenate per-chunk VEP parquets and merge in NMD/AbSplice → annotations.parquet."""


def _vep_chunk_parquets(wildcards):
    return expand(str(OUT / "vep" / "chunk_{chunk}.vep.parquet"), chunk=chunk_ids())


rule merge_vep:
    input:
        parquets=_vep_chunk_parquets,
    output:
        parquet=OUT / "vep.parquet",
    conda:
        "../../envs/parse.yaml"
    run:
        import polars as pl
        pl.concat([pl.scan_parquet(p) for p in input.parquets]).sink_parquet(output.parquet)


rule merge:
    """Final annotation table. NMD (variant×transcript) + AbSplice (variant) folded in if enabled."""
    input:
        vep=OUT / "vep.parquet",
        nmd=(OUT / "nmd.parquet") if NMD else [],
        absplice=(OUT / "absplice" / "absplice.parquet") if ABSPLICE else [],
    output:
        parquet=OUT / "annotations.parquet",
    params:
        with_nmd=NMD,
        with_absplice=ABSPLICE,
    conda:
        "../../envs/parse.yaml"
    script:
        "../scripts/merge_annotations.py"
