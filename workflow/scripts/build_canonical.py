"""Derive the canonical transcript per gene from the gencode GTF.

VEP's `--canonical` is inert in `--gtf` mode: VEP 113 only sets the CANONICAL flag from an
`Ensembl_canonical` GTF tag (Bio::EnsEMBL::VEP::AnnotationSource::File::BaseGXF.pm), and
gencode v34 predates that tag — so VEP's CANONICAL column comes out empty. We reconstruct
it from the GTF, keyed on the versioned transcript id (== VEP `Feature`), using the SAME
rule as the gtex-benchmark consumer (modules/variants/build_transcript_metadata.py) so the
two agree across the repo boundary: canonical = the gene's MANE_Select transcript, else its
best `appris_principal_N`; genes with neither get no canonical.

CANONICAL is a pure property of the transcript (independent of the variant), so `merge`
joins this small table onto the VEP output by `Feature` — no need to re-run VEP.

Output: ref/canonical.parquet (`Feature`, `CANONICAL`="YES", one row per canonical tx).
Invoked by Snakemake (`snakemake.input.gtf`, `snakemake.output.parquet`).
"""
# NB: no `from __future__ import annotations` — Snakemake prepends a preamble to `script:`
# files, which would push a future import off line 1 and raise SyntaxError. Needs py>=3.11.
import polars as pl


def main(gtf: str, out: str) -> None:
    meta = (
        pl.scan_csv(gtf, separator="\t", comment_prefix="#", has_header=False,
                    new_columns=[f"c{i}" for i in range(9)])
        .filter(pl.col("c2") == "transcript")
        .select(
            transcript_id=pl.col("c8").str.extract(r'transcript_id "([^"]+)"', 1),
            gene_id=pl.col("c8").str.extract(r'gene_id "([^"]+)"', 1),
            appris=pl.col("c8").str.extract(r'tag "(appris_[a-z0-9_]+)"', 1),
            mane_select=pl.col("c8").str.contains('tag "MANE_Select"'),
        )
        .with_columns(
            principal=pl.col("appris").str.extract(r"appris_principal_([0-9])", 1).cast(pl.Int8),
        )
        .collect()
    )

    # eligible = MANE_Select or any appris_principal; priority MANE (0) < principal_1..5.
    eligible = pl.col("mane_select") | pl.col("principal").is_not_null()
    priority = pl.when(pl.col("mane_select")).then(0).otherwise(pl.col("principal"))

    canon = (
        meta.filter(eligible)
        .with_columns(_prio=priority)
        .sort("gene_id", "_prio", "transcript_id")            # tie -> smallest transcript id
        .group_by("gene_id", maintain_order=True)
        .agg(Feature=pl.first("transcript_id"))
        .select(pl.col("Feature"), pl.lit("YES").alias("CANONICAL"))
    )
    canon.write_parquet(out)
    print(f"canonical.parquet: {canon.height:,} canonical transcripts "
          f"(of {meta.height:,} total)")


if __name__ == "__main__":
    main(snakemake.input.gtf, snakemake.output.parquet)  # noqa: F821
