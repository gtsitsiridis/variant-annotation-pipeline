"""ENCODE-rE2G enhancer -> gene annotation per variant (hg38, precomputed; no liftOver).

Overlaps each input variant with ENCODE-rE2G enhancer elements and reports the predicted
target gene(s), the rE2G score, and the variant->target-gene-TSS distance. PER-VARIANT only
(no donor/sample/burden aggregation — that is the consumer's job).

Consumes the on-disk E2G export (a normalized star schema):
  - predictions: per-chrom parquet (enhancer_id, target_gene_id[unversioned], target_gene_name,
    target_gene_tss, enhancer_gene_distance, score, model, cell_type_id, chromosome[unprefixed])
  - enhancers.parquet: id (== predictions.enhancer_id), chromosome[unprefixed], start, end, class

Predictions are collapsed to **unique (enhancer, gene) pairs** — the MAX rE2G score across all
biosamples (cell types), with `n_cell_types` recording how many biosamples support the pair.
This dedup happens BEFORE the overlap, so the (expensive) interval range join runs on unique
pairs rather than per-biosample rows. `cell_types` (if set) is an optional pre-filter; null =
all biosamples. Joins predictions -> enhancers on enhancer_id to recover the element interval,
filters by model / score / class, then overlaps the input VCF variants (the VCF is "chr21", the
E2G tables are "21" -> stripped on join; the variant 0-based `start` from the chrom_start_ref_alt
id overlaps the 0-based half-open [enh.start, enh.end)). The unversioned `target_gene_id` is
mapped to the gencode gtf's versioned gene_id. Output is one row per (variant, target_gene)
keyed on variant_id; variants overlapping no enhancer are omitted (absence == not in an enhancer).

Invoked by Snakemake (snakemake.input.vcf; snakemake.params.{predictions,enhancers,gtf,model,
score_threshold,classes,cell_types,memory_limit,threads}; snakemake.output.parquet).
"""

import time
from pathlib import Path

import duckdb
import polars as pl

# Variant x enhancer overlap via interval binning: bucket both sides by floor(pos / BIN_SIZE)
# and hash-equi-join on (chrom, bin) + an exact-overlap filter, instead of a `BETWEEN` range
# join (DuckDB IEJoin) which is ~50x slower on 30M variants. Enhancers are exploded across the
# bins they span; a variant sits in exactly one bin, so no double counting.
BIN_SIZE = 10_000


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


def main(vcf, predictions, enhancers, gtf, out_dir, *, model, score_threshold,
         classes, cell_types, memory_limit="32GB", threads=4):
    t0 = time.time()
    tmp = Path(out_dir).parent / "duckdb_tmp"
    tmp.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    con.execute(f"PRAGMA memory_limit='{memory_limit}'")
    con.execute(f"PRAGMA threads={threads}")
    con.execute("PRAGMA preserve_insertion_order=false")
    con.execute(f"PRAGMA temp_directory='{tmp}'")
    con.register("gene_map", _gene_map(gtf).to_arrow())

    pred_glob = str(Path(predictions) / "*.parquet")
    ct, cl = _sql_in(cell_types), _sql_in(classes)
    # Pre-filter only by model (+ an optional cell-type allowlist). We do NOT keep the biosample
    # dimension: predictions are collapsed to unique (enhancer, gene) pairs with the MAX score
    # across biosamples BEFORE the overlap (so the range join runs on far fewer rows), and the
    # score threshold is applied to that max (HAVING below).
    pred_where = f"model = '{model}'"
    if ct:
        pred_where += f" AND cell_type_id IN {ct}"
    enh_where = f"WHERE class IN {cl}" if cl else ""

    con.execute(f"""
        CREATE TEMP TABLE e2g AS
            WITH preds AS (   -- collapse biosamples: one row per (enhancer, gene), max score
                SELECT enhancer_id, target_gene_id,
                       any_value(target_gene_name) AS target_gene_name,
                       any_value(target_gene_tss) AS target_gene_tss,
                       max(score) AS score,
                       count(DISTINCT cell_type_id) AS n_cell_types
                FROM read_parquet('{pred_glob}')
                WHERE {pred_where}
                GROUP BY enhancer_id, target_gene_id
                HAVING max(score) >= {float(score_threshold)}
            ),
            enh AS (
                SELECT id, chromosome AS chrom, start, "end", class
                FROM read_parquet('{enhancers}') {enh_where}
            ),
            eg AS (   -- enhancer interval + the (biosample-collapsed) gene/score it predicts,
                      -- exploded across the genomic bins it spans (for the equi-join below)
                SELECT e.chrom, e.start, e."end", e.class, p.target_gene_id,
                       p.target_gene_name, p.target_gene_tss, p.score, p.n_cell_types,
                       unnest(range(e.start // {BIN_SIZE},
                                    ((e."end" - 1) // {BIN_SIZE}) + 1)) AS bin
                FROM preds p JOIN enh e ON p.enhancer_id = e.id
            ),
            v AS (    -- input variants: 0-based start parsed from the chrom_start_ref_alt id
                SELECT column2 AS variant_id, replace(column0, 'chr', '') AS chrom,
                       CAST(column1 AS BIGINT) - 1 AS vstart,
                       (CAST(column1 AS BIGINT) - 1) // {BIN_SIZE} AS bin
                FROM read_csv('{vcf}', delim='\t', header=false, comment='#',
                              all_varchar=true, maximum_line_size=10000000)
            ),
            hits AS (   -- hash equi-join on (chrom, bin) then the exact half-open overlap
                SELECT v.variant_id, eg.target_gene_id, eg.target_gene_name,
                       max(eg.score) AS e2g_score,
                       min(abs(v.vstart - eg.target_gene_tss)) AS distance_to_tss,
                       arbitrary(eg.class) AS enhancer_class,
                       max(eg.n_cell_types) AS n_cell_types
                FROM v JOIN eg
                  ON v.chrom = eg.chrom AND v.bin = eg.bin
                  AND v.vstart >= eg.start AND v.vstart < eg."end"
                GROUP BY v.variant_id, eg.target_gene_id, eg.target_gene_name
            )
            SELECT h.variant_id,
                   coalesce(g.gene_id, h.target_gene_id) AS target_gene_id,  -- versioned if mapped
                   h.target_gene_name, h.e2g_score, h.enhancer_class,
                   h.distance_to_tss, h.n_cell_types,
                   -- variant_type (hive partition) from the chrom_start_ref_alt id; e2g is small-only
                   CASE WHEN length(str_split(h.variant_id, '_')[3]) = 1
                         AND length(str_split(h.variant_id, '_')[4]) = 1
                        THEN 'SNV' ELSE 'indel' END AS variant_type
            FROM hits h
            LEFT JOIN gene_map g ON h.target_gene_id = g.gene_id_unversioned
    """)
    for vt in ("SNV", "indel"):
        d = Path(out_dir) / f"variant_type={vt}"
        d.mkdir(parents=True, exist_ok=True)
        con.execute(f"COPY (SELECT * EXCLUDE (variant_type) FROM e2g WHERE variant_type = '{vt}') "
                    f"TO '{d / 'data.parquet'}' (FORMAT parquet)")
    s = con.execute(
        f"SELECT count(*) n, count(DISTINCT variant_id) v, count(DISTINCT target_gene_id) g "
        f"FROM read_parquet('{Path(out_dir)}/**/*.parquet', hive_partitioning=true)"
    ).fetchone()
    print(f"wrote e2g.parquet/variant_type={{SNV,indel}} in {time.time()-t0:.1f}s")
    print(f"  rows={s[0]:,} variants_in_enhancer={s[1]:,} target_genes={s[2]:,}")


if __name__ == "__main__":
    smk = snakemake  # noqa: F821
    p = smk.params
    main(
        smk.input.vcf, p.predictions, p.enhancers, p.gtf, p.out_dir,
        model=p.model, score_threshold=p.score_threshold, classes=p.classes,
        cell_types=p.cell_types, memory_limit=p.memory_limit, threads=int(p.threads),
    )
