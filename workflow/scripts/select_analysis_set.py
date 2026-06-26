"""Pick the deep-tier variant subset from the broad fastVEP (small-track) annotation.

The deep tools (VEP plugins, NMD, AbSplice) are expensive and only meaningful near genes, so
we keep a variant if ANY of its transcript rows is within the deep window: distance <= window,
or distance IS NULL (the variant overlaps the transcript), or it has a non-MODIFIER
consequence. Writes one variant_id per line — the next rule subsets the input VCF on these ids.

Invoked by Snakemake (snakemake.input.basic, snakemake.output.ids, snakemake.params.window).
"""

import polars as pl


def main(basic: str, out_ids: str, window: int) -> None:
    ids = (
        pl.scan_parquet(basic)
        .filter(
            pl.col("distance").is_null()
            | (pl.col("distance") <= window)
            | (pl.col("impact") != "MODIFIER")
        )
        .select("variant_id")
        .unique()
        .collect()
    )
    ids.write_csv(out_ids, include_header=False)
    print(f"analysis set: {ids.height:,} variants (deep window = {window} bp)")


if __name__ == "__main__":
    main(snakemake.input.basic, snakemake.output.ids, int(snakemake.params.window))  # noqa: F821
