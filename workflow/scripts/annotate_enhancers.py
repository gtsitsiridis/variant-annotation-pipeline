"""ENCODE-rE2G enhancer -> gene annotation per variant (hg38, precomputed; no liftOver).

Overlaps each input variant with ENCODE-rE2G enhancer elements and reports the predicted
target gene(s), the rE2G score, and the variant->target-gene-TSS distance. PER-VARIANT only
(no donor/sample/burden aggregation — that is the consumer's job).

Consumes the on-disk E2G export (a normalized star schema):
  - predictions: per-chrom parquet (enhancer_id, target_gene_id[unversioned], target_gene_name,
    target_gene_tss, enhancer_gene_distance, score, model, cell_type_id, chromosome[unprefixed])
  - enhancers.parquet: id (== predictions.enhancer_id), chromosome[unprefixed], start, end, class

Joins predictions -> enhancers on enhancer_id to recover the element interval, filters by
model / score / class / cell type, then overlaps the input VCF variants (the VCF is "chr21",
the E2G tables are "21" -> stripped on join; the variant 0-based `start` from the
chrom_start_ref_alt id overlaps the 0-based half-open [enh.start, enh.end)). The unversioned
`target_gene_id` is mapped to the gencode gtf's versioned gene_id. Output is one row per
(variant, target_gene) keyed on variant_id, with the max score across the selected cell types;
variants overlapping no enhancer are omitted (absence == not in an enhancer).

Invoked by Snakemake (snakemake.input.vcf; snakemake.params.{predictions,enhancers,gtf,model,
score_threshold,classes,cell_types,memory_limit,threads}; snakemake.output.parquet).
"""
from __future__ import annotations

import time
from pathlib import Path

import duckdb
import polars as pl


def _gene_map(gtf: str) -> pl.DataFrame:
    """unversioned ENSG -> versioned gene_id, from the gencode GTF gene lines."""
    return (
        pl.scan_csv(gtf, separator="\t", comment_prefix="#", has_header=False,
                    new_columns=[f"c{i}" for i in range(9)])
        .filter(pl.col("c2") == "gene")
        .select(gene_id=pl.col("c8").str.extract(r'gene_id "([^"]+)"', 1))
        .with_columns(gene_id_unversioned=pl.col("gene_id").str.replace(r"\..*$", ""))
        .unique(subset="gene_id_unversioned")
        .collect()
    )


def _sql_in(vals) -> str | None:
    """A SQL IN-list from a python list, or None to skip the filter."""
    if not vals:
        return None
    return "(" + ",".join("'" + str(v).replace("'", "''") + "'" for v in vals) + ")"


def main(vcf, predictions, enhancers, gtf, out, *, model, score_threshold,
         classes, cell_types, memory_limit="32GB", threads=4):
    t0 = time.time()
    tmp = Path(out).parent / "duckdb_tmp"
    tmp.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    con.execute(f"PRAGMA memory_limit='{memory_limit}'")
    con.execute(f"PRAGMA threads={threads}")
    con.execute("PRAGMA preserve_insertion_order=false")
    con.execute(f"PRAGMA temp_directory='{tmp}'")
    con.register("gene_map", _gene_map(gtf).to_arrow())

    pred_glob = str(Path(predictions) / "*.parquet")
    ct, cl = _sql_in(cell_types), _sql_in(classes)
    pred_where = f"model = '{model}' AND score >= {float(score_threshold)}"
    if ct:
        pred_where += f" AND cell_type_id IN {ct}"
    enh_where = f"WHERE class IN {cl}" if cl else ""

    con.execute(f"""
        COPY (
            WITH preds AS (
                SELECT enhancer_id, target_gene_id, target_gene_name, target_gene_tss,
                       score, cell_type_id
                FROM read_parquet('{pred_glob}')
                WHERE {pred_where}
            ),
            enh AS (
                SELECT id, chromosome AS chrom, start, "end", class
                FROM read_parquet('{enhancers}') {enh_where}
            ),
            eg AS (   -- enhancer interval + the gene/score it predicts
                SELECT e.chrom, e.start, e."end", e.class, p.target_gene_id,
                       p.target_gene_name, p.target_gene_tss, p.score, p.cell_type_id
                FROM preds p JOIN enh e ON p.enhancer_id = e.id
            ),
            v AS (    -- input variants: 0-based start parsed from the chrom_start_ref_alt id
                SELECT column2 AS variant_id, replace(column0, 'chr', '') AS chrom,
                       CAST(column1 AS BIGINT) - 1 AS vstart
                FROM read_csv('{vcf}', delim='\t', header=false, comment='#',
                              all_varchar=true, maximum_line_size=10000000)
            ),
            hits AS (
                SELECT v.variant_id, eg.target_gene_id, eg.target_gene_name,
                       max(eg.score) AS e2g_score,
                       min(abs(v.vstart - eg.target_gene_tss)) AS distance_to_tss,
                       arbitrary(eg.class) AS enhancer_class,
                       count(DISTINCT eg.cell_type_id) AS n_cell_types
                FROM v JOIN eg
                  ON v.chrom = eg.chrom AND v.vstart >= eg.start AND v.vstart < eg."end"
                GROUP BY v.variant_id, eg.target_gene_id, eg.target_gene_name
            )
            SELECT h.variant_id,
                   coalesce(g.gene_id, h.target_gene_id) AS target_gene_id,  -- versioned if mapped
                   h.target_gene_name, h.e2g_score, h.enhancer_class,
                   h.distance_to_tss, h.n_cell_types
            FROM hits h
            LEFT JOIN gene_map g ON h.target_gene_id = g.gene_id_unversioned
        ) TO '{out}' (FORMAT parquet)
    """)
    s = con.execute(f"""
        SELECT count(*) n, count(DISTINCT variant_id) v, count(DISTINCT target_gene_id) g
        FROM read_parquet('{out}')
    """).fetchone()
    print(f"wrote {Path(out).name} in {time.time()-t0:.1f}s")
    print(f"  rows={s[0]:,} variants_in_enhancer={s[1]:,} target_genes={s[2]:,}")


if __name__ == "__main__":
    smk = snakemake  # noqa: F821
    p = smk.params
    main(
        smk.input.vcf, p.predictions, p.enhancers, p.gtf, smk.output.parquet,
        model=p.model, score_threshold=p.score_threshold, classes=p.classes,
        cell_types=p.cell_types, memory_limit=p.memory_limit, threads=int(p.threads),
    )
