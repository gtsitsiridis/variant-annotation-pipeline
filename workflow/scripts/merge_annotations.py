"""Merge the concatenated VEP annotations with NMD-Scanner and AbSplice → final table.

VEP and NMD are transcript-level (one row per variant×transcript); AbSplice is per variant.
- NMD joins on the variant key **+ transcript** (`Feature` == NMD `transcript_id`).
- AbSplice joins on the variant key only (broadcast across the variant's transcripts).

Invoked by Snakemake (`snakemake.input.vep`, optional `.nmd`/`.absplice`,
`snakemake.params.with_nmd`/`with_absplice`, `snakemake.output.parquet`).
"""
from __future__ import annotations

import polars as pl

KEY = ["chrom", "start", "ref", "alt"]


def _one(path) -> str | None:
    """Snakemake passes [] for a disabled optional input; normalize to a path or None."""
    if isinstance(path, list):
        return path[0] if path else None
    return path or None


def main(vep: str, out: str, nmd: str | None, absplice: str | None) -> None:
    lf = pl.scan_parquet(vep)
    if nmd:
        lf = lf.join(pl.scan_parquet(nmd), on=KEY + ["Feature"], how="left")
    if absplice:
        lf = lf.join(pl.scan_parquet(absplice), on=KEY, how="left")
    lf.sink_parquet(out)


if __name__ == "__main__":
    smk = snakemake  # noqa: F821
    main(
        smk.input.vep,
        smk.output.parquet,
        _one(smk.input.nmd) if smk.params.with_nmd else None,
        _one(smk.input.absplice) if smk.params.with_absplice else None,
    )
