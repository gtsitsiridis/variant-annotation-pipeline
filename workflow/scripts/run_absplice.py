"""Aggregate an AbSplice-DNA result table to a per-variant annotation parquet.

AbSplice (gagneurlab/absplice) is run separately (its own snakemake + per-tissue
SpliceMaps); this consumes its output and produces a tidy per-variant table:
    chrom, start, ref, alt, AbSplice_DNA            # max across the requested tissues
    [+ AbSplice_DNA_<tissue> ...]                   # if --per-tissue true

The AbSplice output is expected to have a variant id, a tissue, and an `AbSplice_DNA`
score (column names vary by AbSplice version; we detect common ones). Variant ids like
`chr17:41197805:T>C` or `17:41197805:T:C` are parsed; AbSplice pos is 1-based (VCF) →
converted to our 0-based `start` (= pos - 1) to match the variant key.
"""
from __future__ import annotations

import argparse

import polars as pl

# Candidate column names across AbSplice versions.
VARIANT_COLS = ["variant", "variant_id", "var", "variant_name"]
TISSUE_COLS = ["tissue", "Tissue"]
SCORE_COLS = ["AbSplice_DNA", "absplice_dna", "AbSplice_DNA_score"]


def _pick(cols: list[str], candidates: list[str], what: str) -> str:
    for c in candidates:
        if c in cols:
            return c
    raise SystemExit(f"AbSplice result has no {what} column (looked for {candidates}); got {cols}")


def parse_variant(expr: pl.Expr) -> dict[str, pl.Expr]:
    """`chr17:41197805:T>C` / `17:41197805:T:C` → chrom, start(0-based), ref, alt."""
    norm = expr.str.replace(">", ":", literal=True)  # unify the ref/alt separator
    parts = norm.str.split(":")
    chrom = parts.list.get(0)
    chrom = pl.when(chrom.str.starts_with("chr")).then(chrom).otherwise("chr" + chrom)
    return {
        "chrom": chrom,
        "start": parts.list.get(1).cast(pl.Int64, strict=False) - 1,  # 1-based → 0-based
        "ref": parts.list.get(2),
        "alt": parts.list.get(3),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--absplice-result", required=True, help="AbSplice output (.parquet/.csv/.tsv)")
    ap.add_argument("--tissues", default="", help="comma-separated; empty = all tissues")
    ap.add_argument("--per-tissue", default="false")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    path = args.absplice_result
    lf = pl.scan_parquet(path) if path.endswith((".parquet", ".pq")) else pl.scan_csv(
        path, separator="\t" if path.endswith(".tsv") else ","
    )
    cols = lf.collect_schema().names()
    vcol = _pick(cols, VARIANT_COLS, "variant")
    tcol = _pick(cols, TISSUE_COLS, "tissue")
    scol = _pick(cols, SCORE_COLS, "AbSplice_DNA")

    df = lf.with_columns(**parse_variant(pl.col(vcol))).rename({scol: "AbSplice_DNA", tcol: "tissue"})
    tissues = [t for t in args.tissues.split(",") if t]
    if tissues:
        df = df.filter(pl.col("tissue").is_in(tissues))

    key = ["chrom", "start", "ref", "alt"]
    agg = df.group_by(key).agg(pl.col("AbSplice_DNA").max())

    if args.per_tissue.lower() == "true":
        wide = (
            df.collect()
            .pivot(values="AbSplice_DNA", index=key, on="tissue", aggregate_function="max")
            .rename(lambda c: f"AbSplice_DNA_{c}" if c not in key else c)
        )
        agg = agg.collect().join(wide, on=key, how="left").lazy()

    agg.sink_parquet(args.out)


if __name__ == "__main__":
    main()
