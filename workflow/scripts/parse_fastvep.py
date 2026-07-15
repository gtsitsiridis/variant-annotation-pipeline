"""fastVEP VCF/CSQ output -> tidy transcript-level parquet, hive-partitioned by variant_type.

Explodes INFO/CSQ to one row per variant x transcript (all transcripts kept), keyed on the
VCF ID column (`variant_id` = chrom_start_ref_alt for small variants, sv_id for SVs — fastVEP
round-trips whatever ID it was given). Joins transcript_metadata (symbol/biotype/canonical/tsl)
on the transcript id. Frequency columns are NOT joined here (dataset-specific; the consumer
joins them back on `variant_id`).

Output layout: <out_dir>/variant_type=<vtype>/<track>.parquet, where
  track "small" -> variant_type in {SNV, indel} (split by ref/alt length),
  track "sv"    -> variant_type = SV.
So the union `<out_dir>/**/*.parquet` (hive_partitioning=true) is the single `fastvep` table
partitioned by variant_type — the same layout every tool in this pipeline writes.

CSQ field order (1-based, see ##INFO CSQ header): 2 Consequence, 3 IMPACT, 5 Gene, 7 Feature,
9 EXON, 10 INTRON, 11 HGVSc, 12 HGVSp, 13 cDNA, 14 CDS, 15 Protein, 16 Amino_acids, 17 Codons,
21 DISTANCE, 22 STRAND, 33 ENSP.  Intergenic rows have Gene='-' and are dropped.

DuckDB (not polars) because the small track explodes to ~110M rows.

Invoked by Snakemake (snakemake.input.vcf/.metadata, snakemake.params.out_dir/.track/
.memory_limit/.threads).
"""

import time
from pathlib import Path

import duckdb

# which variant_type partition(s) a track produces, and the SQL predicate selecting each
_TRACK_PARTITIONS = {
    "small": {"SNV": "length(ref) = 1 AND length(alt) = 1",
              "indel": "NOT (length(ref) = 1 AND length(alt) = 1)"},
    "sv": {"SV": "TRUE"},
}


def parse(vcf: str, meta: str, out_dir: str, track: str,
          memory_limit: str = "32GB", threads: int = 4) -> None:
    t0 = time.time()
    partitions = _TRACK_PARTITIONS[track]
    out_dir = Path(out_dir)
    tmp = out_dir / "duckdb_tmp"
    tmp.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    con.execute(f"PRAGMA memory_limit='{memory_limit}'")
    con.execute(f"PRAGMA threads={threads}")
    con.execute("PRAGMA preserve_insertion_order=false")
    con.execute(f"PRAGMA temp_directory='{tmp}'")
    read_vcf = (
        f"read_csv('{vcf}', delim='\t', header=false, comment='#', "
        f"all_varchar=true, maximum_line_size=10000000)"
    )

    # materialize the transcript-level annotation once, then write one file per variant_type
    con.execute(f"""
        CREATE TEMP TABLE final AS
        WITH recs AS (
            SELECT column2 AS variant_id, column0 AS chrom,
                   CAST(column1 AS BIGINT) - 1 AS "start", column3 AS ref, column4 AS alt,
                   regexp_extract(column7, 'CSQ=([^;]*)', 1) AS csq
            FROM {read_vcf}
        ),
        ex AS (
            SELECT variant_id, chrom, "start", ref, alt,
                   string_split(unnest(string_split(csq, ',')), '|') AS p
            FROM recs WHERE csq <> ''
        ),
        ann AS (
            SELECT variant_id, chrom, "start", ref, alt,
                   p[5] AS gene, p[7] AS transcript, p[2] AS consequence, p[3] AS impact,
                   CASE WHEN p[21]='' THEN NULL ELSE CAST(p[21] AS INTEGER) END AS distance,
                   NULLIF(p[9],'')  AS exon,  NULLIF(p[10],'') AS intron,
                   NULLIF(p[11],'') AS hgvsc, NULLIF(p[12],'') AS hgvsp,
                   NULLIF(p[13],'') AS cdna_position, NULLIF(p[14],'') AS cds_position,
                   NULLIF(p[15],'') AS protein_position, NULLIF(p[16],'') AS amino_acids,
                   NULLIF(p[17],'') AS codons,
                   CASE WHEN p[22]='' THEN NULL ELSE CAST(p[22] AS INTEGER) END AS strand,
                   NULLIF(p[33],'') AS ensp
            FROM ex WHERE p[5] <> '-'   -- drop intergenic (Gene='-')
        )
        SELECT a.variant_id, a.chrom, a."start", a.ref, a.alt,
               a.gene, m.symbol, a.transcript, m.biotype, m.canonical, m.tsl,
               a.consequence, a.impact, a.distance, a.exon, a.intron,
               a.hgvsc, a.hgvsp, a.ensp,
               a.cdna_position, a.cds_position, a.protein_position,
               a.amino_acids, a.codons, a.strand
        FROM ann a
        LEFT JOIN read_parquet('{meta}') m ON a.transcript = m.transcript_id
    """)

    for vtype, where in partitions.items():
        d = out_dir / f"variant_type={vtype}"
        d.mkdir(parents=True, exist_ok=True)
        con.execute(f"COPY (SELECT * FROM final WHERE {where}) "
                    f"TO '{d / (track + '.parquet')}' (FORMAT parquet)")

    s = con.execute("""
        SELECT count(*) n_rows, count(DISTINCT variant_id) n_variants,
               count(DISTINCT transcript) n_tx, count(DISTINCT gene) n_genes,
               count(*) FILTER (WHERE canonical) n_canon,
               count(*) FILTER (WHERE hgvsc IS NOT NULL) n_hgvsc,
               count(*) FILTER (WHERE symbol IS NULL) n_no_meta
        FROM final
    """).fetchone()
    print(f"wrote fastvep/{{{','.join(partitions)}}} for track={track} in {time.time()-t0:.1f}s")
    print(f"  rows={s[0]:,} variants={s[1]:,} transcripts={s[2]:,} genes={s[3]:,} "
          f"canonical={s[4]:,} hgvsc={s[5]:,} no_meta={s[6]:,}")


if __name__ == "__main__":
    parse(
        snakemake.input.vcf, snakemake.input.metadata,                    # noqa: F821
        snakemake.params.out_dir, snakemake.params.track,                 # noqa: F821
        snakemake.params.memory_limit, int(snakemake.params.threads),     # noqa: F821
    )
