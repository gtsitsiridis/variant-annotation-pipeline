"""gencode GTF -> transcript metadata for joining onto fastVEP annotations.

fastVEP run on the gencode GFF3 does NOT populate SYMBOL/BIOTYPE/CANONICAL/MANE/APPRIS, so
we extract them from the gencode GTF, keyed on the versioned transcript id (ENST...N), which
matches fastVEP's `Feature` column. `canonical` is derived (gencode v34 has no
Ensembl_canonical tag): MANE_Select where present, else the gene's best appris_principal_N;
genes with neither get no canonical.

Invoked by Snakemake (snakemake.input.gtf, snakemake.output.parquet).
"""

import polars as pl


def build(gtf: str, out: str) -> None:
    # read GTF, keep transcript rows, pull the attribute column
    lf = (
        pl.scan_csv(gtf, separator="\t", comment_prefix="#", has_header=False,
                    new_columns=[f"c{i}" for i in range(9)])
        .filter(pl.col("c2") == "transcript")
        .select(attr=pl.col("c8"))
    )
    g = lambda pat: pl.col("attr").str.extract(pat, 1)  # noqa: E731
    meta = lf.select(
        transcript_id=g(r'transcript_id "([^"]+)"'),
        gene_id=g(r'gene_id "([^"]+)"'),
        symbol=g(r'gene_name "([^"]+)"'),
        gene_type=g(r'gene_type "([^"]+)"'),
        biotype=g(r'transcript_type "([^"]+)"'),
        tsl=g(r'transcript_support_level "([^"]+)"'),
        appris=g(r'tag "(appris_[a-z0-9_]+)"'),
        mane_select=pl.col("attr").str.contains('tag "MANE_Select"'),
        basic=pl.col("attr").str.contains('tag "basic"'),
        ccds=pl.col("attr").str.contains('tag "CCDS"'),
    ).collect()

    # principal rank for the canonical fallback
    meta = meta.with_columns(
        principal=pl.col("appris").str.extract(r"appris_principal_([0-9])", 1).cast(pl.Int8),
    )
    # priority among "eligible" transcripts: MANE first (0), else principal_1..5
    eligible = pl.col("mane_select") | pl.col("principal").is_not_null()
    priority = pl.when(pl.col("mane_select")).then(0).otherwise(pl.col("principal"))

    # one canonical transcript per gene = lowest priority among eligible (tie -> first id)
    canon = (
        meta.filter(eligible)
        .with_columns(_prio=priority)
        .sort("gene_id", "_prio", "transcript_id")
        .group_by("gene_id", maintain_order=True)
        .agg(canonical_tx=pl.first("transcript_id"))
    )
    meta = (
        meta.join(canon, on="gene_id", how="left")
        .with_columns(
            canonical=(pl.col("transcript_id") == pl.col("canonical_tx")).fill_null(False)
        )
        .drop("canonical_tx")
        .with_columns(tsl=pl.when(pl.col("tsl") == "NA").then(None)
                      .otherwise(pl.col("tsl")).cast(pl.Int8, strict=False))
    )

    meta.write_parquet(out)
    print(f"transcript_metadata.parquet: {meta.height:,} transcripts")
    print(f"  canonical: {meta['canonical'].sum():,}  | MANE_Select: {meta['mane_select'].sum():,}")
    print(f"  symbol non-null: {meta['symbol'].is_not_null().sum():,}")


if __name__ == "__main__":
    build(snakemake.input.gtf, snakemake.output.parquet)  # noqa: F821
