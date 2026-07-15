"""Parse a VEP `--tab` output (one chromosome) into a tidy parquet.

VEP tab has a `##`-prefixed header block then a single `#Uploaded_variation`-prefixed
column header. `Uploaded_variation` is our join key `chrom_start_ref_alt`. Output is
transcript-level (one row per variant×transcript) with the numeric plugin scores parsed
out and the consequence terms one-hot encoded. Column selection is defensive — VEP/plugin
versions differ, so we keep whatever is present.

Invoked by Snakemake (`snakemake.input.tab`, `snakemake.output.parquet`).
"""
# NB: no `from __future__ import annotations` — Snakemake prepends a preamble to `script:`
# files, which would push a future import off line 1 and raise SyntaxError. Needs py>=3.11.
import polars as pl

# Consequence terms one-hot encoded (the DeepRVAT set; HIGH/MODERATE coding + splicing).
CONSEQUENCES = [
    "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant",
    "stop_gained", "stop_lost", "start_lost", "frameshift_variant",
    "inframe_insertion", "inframe_deletion", "missense_variant",
    "protein_altering_variant", "synonymous_variant",
    "5_prime_UTR_variant", "3_prime_UTR_variant",
    "upstream_gene_variant", "downstream_gene_variant", "intron_variant",
]
SPLICEAI_DS = ["SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL",
               "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL"]


def _score_in_parens(col: str) -> pl.Expr:
    """VEP SIFT/PolyPhen come as 'deleterious(0.03)' → extract the float."""
    return pl.col(col).str.extract(r"\(([0-9.eE+-]+)\)", 1).cast(pl.Float64, strict=False)


def main(tab: str, out: str) -> None:
    lf = pl.scan_csv(
        tab, separator="\t", comment_prefix="##", null_values=["-", ""],
        infer_schema_length=10000, truncate_ragged_lines=True,
    )
    cols = lf.collect_schema().names()
    # First column is '#Uploaded_variation'; normalize.
    upload = cols[0]
    lf = lf.rename({upload: "Uploaded_variation"})

    # Variant key from the ID we wrote upstream: chrom_start_ref_alt.
    key = pl.col("Uploaded_variation").str.splitn("_", 4)
    lf = lf.with_columns(
        chrom=key.struct[0],
        start=key.struct[1].cast(pl.Int64, strict=False),
        ref=key.struct[2],
        alt=key.struct[3],
    )

    have = set(cols)
    out_cols = [
        pl.col("Uploaded_variation").alias("variant_id"),   # join key (chrom_start_ref_alt)
        pl.col("chrom"), pl.col("start"), pl.col("ref"), pl.col("alt"),
    ]
    # Transcript-level categorical fields (keep if present).
    for c in ["Gene", "Feature", "BIOTYPE", "CANONICAL", "IMPACT", "SYMBOL", "Consequence"]:
        if c in have:
            out_cols.append(pl.col(c))
    # SIFT / PolyPhen → numeric score.
    if "SIFT" in have:
        out_cols.append(_score_in_parens("SIFT").alias("sift_score"))
    if "PolyPhen" in have:
        out_cols.append(_score_in_parens("PolyPhen").alias("polyphen_score"))
    # Plugin numeric scores (names per VEP plugin output).
    for src, dst in [
        ("CADD_RAW", "CADD_raw"), ("CADD_PHRED", "CADD_PHRED"),
        ("PrimateAI", "PrimateAI_score"), ("am_pathogenicity", "alphamissense"),
        ("Condel", "Condel"), ("gnomADg_AF", "gnomADg_AF"), ("gnomADe_AF", "gnomADe_AF"),
    ]:
        if src in have:
            out_cols.append(pl.col(src).cast(pl.Float64, strict=False).alias(dst))
    # SpliceAI: max delta across DS_AG/AL/DG/DL. VEP's SpliceAI plugin emits one combined
    # column `SpliceAI_pred` = SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL (the 4
    # delta scores are fields 1-4); older setups emit four separate SpliceAI_pred_DS_* cols.
    if "SpliceAI_pred" in have:
        ds = pl.col("SpliceAI_pred").str.split("|")
        out_cols.append(
            pl.max_horizontal([
                ds.list.get(i, null_on_oob=True).cast(pl.Float64, strict=False)
                for i in (1, 2, 3, 4)
            ]).alias("SpliceAI_delta_score")
        )
    else:
        ds_present = [c for c in SPLICEAI_DS if c in have]
        if ds_present:
            out_cols.append(
                pl.max_horizontal([pl.col(c).cast(pl.Float64, strict=False) for c in ds_present])
                .alias("SpliceAI_delta_score")
            )
    # LOFTEE.
    for c in ["LoF", "LoF_filter", "LoF_flags"]:
        if c in have:
            out_cols.append(pl.col(c))

    df = lf.select(out_cols).with_columns(   # variant_type (for the hive partition) from ref/alt
        variant_type=pl.when((pl.col("ref").str.len_chars() == 1) & (pl.col("alt").str.len_chars() == 1))
                       .then(pl.lit("SNV")).otherwise(pl.lit("indel"))
    )
    # One-hot consequence terms (Consequence is '&'-joined).
    if "Consequence" in have:
        for term in CONSEQUENCES:
            df = df.with_columns(
                pl.col("Consequence").str.contains(term, literal=True).cast(pl.Int8)
                .alias(f"Consequence_{term}")
            )

    df.sink_parquet(out)


if __name__ == "__main__":
    main(snakemake.input.tab, snakemake.output.parquet)  # noqa: F821
