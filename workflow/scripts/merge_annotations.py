"""Merge the concatenated VEP annotations with AbSplice into the final table.

VEP is transcript-level (one row per variant×transcript); AbSplice is per variant
(aggregated across tissues by `run_absplice.py`). We left-join AbSplice onto the VEP rows
on the variant key so every transcript row carries the variant's AbSplice score.

Invoked by Snakemake (`snakemake.input.vep`, optional `snakemake.input.absplice`,
`snakemake.params.with_absplice`, `snakemake.output.parquet`).
"""
from __future__ import annotations

import polars as pl

KEY = ["chrom", "start", "ref", "alt"]


def main(vep: str, out: str, absplice: str | None) -> None:
    lf = pl.scan_parquet(vep)
    if absplice:
        lf = lf.join(pl.scan_parquet(absplice), on=KEY, how="left")
    lf.sink_parquet(out)


if __name__ == "__main__":
    smk = snakemake  # noqa: F821
    abs_path = smk.input.absplice if smk.params.with_absplice else None
    if isinstance(abs_path, list):
        abs_path = abs_path[0] if abs_path else None
    main(smk.input.vep, smk.output.parquet, abs_path)
