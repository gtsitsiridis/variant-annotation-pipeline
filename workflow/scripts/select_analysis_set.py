"""Filter the fastVEP transcript table -> the variant ids that pass the pipeline filters.

fastVEP is the funnel: a small variant is kept for the `additional` tools if ANY of its transcript
rows is within the cis window (distance <= window, or distance IS NULL = the variant overlaps the
transcript) AND — when the corresponding flag is on — the transcript's gene is protein_coding AND
the transcript is canonical. Writes one variant_id per line; the next rule subsets input_vcf on
these ids into fastvep/small.vcf.gz (the VCF every additional tool re-annotates).

Invoked by Snakemake (snakemake.input.fastvep; snakemake.output.ids;
snakemake.params.window/.protein_coding_only/.canonical_only).
"""

import polars as pl


def main(fastvep_parquet: str, out_ids: str, window: int,
         protein_coding_only: bool, canonical_only: bool) -> None:
    lf = (
        pl.scan_parquet(fastvep_parquet)
        .filter(pl.col("variant_type").is_in(["SNV", "indel"]))   # small track only
        .filter(pl.col("distance").is_null() | (pl.col("distance") <= window))
    )
    if protein_coding_only:
        lf = lf.filter(pl.col("gene_type") == "protein_coding")
    if canonical_only:
        lf = lf.filter(pl.col("canonical"))
    ids = lf.select("variant_id").unique().collect()
    ids.write_csv(out_ids, include_header=False)
    print(f"fastVEP-filtered analysis set: {ids.height:,} variants "
          f"(window={window}, protein_coding_only={protein_coding_only}, "
          f"canonical_only={canonical_only})")


if __name__ == "__main__":
    smk = snakemake  # noqa: F821
    main(smk.input.fastvep, smk.output.ids, int(smk.params.window),
         bool(smk.params.protein_coding_only), bool(smk.params.canonical_only))
