"""Per-(variant, gene) AbSplice2 scores from the precomputed gagneurlab/absplice result.

We do NOT run AbSplice (its own heavy workflow needs per-tissue SpliceMaps + SpliceAI +
GPU); we consume the precomputed result deeprvat uses: one parquet *dataset per gene*
(≈19k `ENSG…_max_preds.parquet/` dirs), keyed on (chrom, start[0-based], ref, alt,
gene_id) — **SNV-only** — carrying `AbSplice_DNA_max` and `AbSplice2_max` (max across the
49 GTEx tissues) plus per-tissue columns.

We scan it all, project the two max scores, and inner-join to our variant×gene keys (from
vep.parquet; the gencode version is stripped from VEP's `Gene` to match AbSplice's
`gene_id`). The merge step then attaches these on (chrom,start,ref,alt,gene); indels —
absent from the SNV-only source — get null.

Output: chrom, start, ref, alt, gene, AbSplice_DNA_max, AbSplice2_max
Invoked by Snakemake (`snakemake.input.absplice_dir`, `.vep`, `snakemake.output.parquet`).
"""
# NB: no `from __future__ import annotations` — Snakemake prepends a `script:` preamble.
import os

import polars as pl

KEY = ["chrom", "start", "ref", "alt", "gene"]


def main(absplice_dir: str, vep: str, out: str) -> None:
    # our variant×gene keys (strip the gencode version: ENSG….N -> ENSG…)
    keys = (
        pl.scan_parquet(vep)
        .select("chrom", "start", "ref", "alt",
                pl.col("Gene").str.replace(r"\..*$", "").alias("gene"))
        .unique()
    )
    # precomputed AbSplice2 (cat columns -> str); start is already 0-based like our key.
    # Gene files vary in int width (start Int32 vs Int64) -> upcast on scan.
    absp = (
        pl.scan_parquet(os.path.join(absplice_dir, "**", "*.parquet"), hive_partitioning=False,
                        cast_options=pl.ScanCastOptions(integer_cast="upcast"))
        .select(pl.col("chrom").cast(pl.String), pl.col("start").cast(pl.Int64),
                pl.col("ref").cast(pl.String), pl.col("alt").cast(pl.String),
                pl.col("gene_id").cast(pl.String).alias("gene"),
                "AbSplice_DNA_max", "AbSplice2_max")
    )
    keys.join(absp, on=KEY, how="inner").sink_parquet(out)


if __name__ == "__main__":
    smk = snakemake  # noqa: F821
    main(smk.input.absplice_dir, smk.input.vep, smk.output.parquet)
