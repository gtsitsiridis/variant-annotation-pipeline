"""fastVEP CSQ VCF(s) -> ONE flat transcript-level parquet (parts/fastvep.parquet).

Explodes INFO/CSQ to one row per variant × transcript (all transcripts kept), keyed on the
VCF ID column (`variant_id` = chrom_start_ref_alt for small variants, sv_id for SVs — fastVEP
round-trips whatever ID it was given). Both tracks (small + optional sv) are UNIONed into one
table with a `variant_type` column (SNV/indel from ref/alt length; SV for the sv track). Joins
transcript_metadata (gene_type/symbol/biotype/canonical/tsl) on the transcript id — the
`gene_type` + `canonical` columns drive the downstream filter/partitioning in combine.

No hive partitioning here (that is combine's job): this is a single flat parquet FILE that both
`fastvep_filtered_vcf` (to pick the analysis set) and `combine_annotations` read.

CSQ field order (1-based, see ##INFO CSQ header): 2 Consequence, 3 IMPACT, 5 Gene, 7 Feature,
9 EXON, 10 INTRON, 11 HGVSc, 12 HGVSp, 13 cDNA, 14 CDS, 15 Protein, 16 Amino_acids, 17 Codons,
21 DISTANCE, 22 STRAND, 33 ENSP.  Intergenic rows have Gene='-' and are dropped.

DuckDB (not polars) because the small track explodes to ~110M rows.

Invoked by Snakemake (snakemake.input.small [+ .sv] / .metadata; snakemake.output.parquet;
snakemake.params.memory_limit/.threads).
"""

import time
from pathlib import Path

import duckdb


def _track_sql(vcf: str, track: str) -> str:
    """Exploded transcript rows + a variant_type expression for one CSQ VCF track."""
    read_vcf = (
        f"read_csv('{vcf}', delim='\t', header=false, comment='#', "
        f"all_varchar=true, maximum_line_size=10000000)"
    )
    vtype = ("'SV'" if track == "sv"
             else "CASE WHEN length(ref) = 1 AND length(alt) = 1 THEN 'SNV' ELSE 'indel' END")
    return f"""
        SELECT variant_id, chrom, "start", ref, alt, {vtype} AS variant_type,
               p[5] AS gene, p[7] AS transcript, p[2] AS consequence, p[3] AS impact,
               CASE WHEN p[21]='' THEN NULL ELSE CAST(p[21] AS INTEGER) END AS distance,
               NULLIF(p[9],'')  AS exon,  NULLIF(p[10],'') AS intron,
               NULLIF(p[11],'') AS hgvsc, NULLIF(p[12],'') AS hgvsp,
               NULLIF(p[13],'') AS cdna_position, NULLIF(p[14],'') AS cds_position,
               NULLIF(p[15],'') AS protein_position, NULLIF(p[16],'') AS amino_acids,
               NULLIF(p[17],'') AS codons,
               CASE WHEN p[22]='' THEN NULL ELSE CAST(p[22] AS INTEGER) END AS strand,
               NULLIF(p[33],'') AS ensp
        FROM (
            SELECT variant_id, chrom, "start", ref, alt,
                   string_split(unnest(string_split(csq, ',')), '|') AS p
            FROM (
                SELECT column2 AS variant_id, column0 AS chrom,
                       CAST(column1 AS BIGINT) - 1 AS "start", column3 AS ref, column4 AS alt,
                       regexp_extract(column7, 'CSQ=([^;]*)', 1) AS csq
                FROM {read_vcf}
            ) WHERE csq <> ''
        ) WHERE p[5] <> '-'   -- drop intergenic (Gene='-')
    """


def parse(small_vcf: str, sv_vcf, meta: str, out_parquet: str,
          memory_limit: str = "32GB", threads: int = 4) -> None:
    t0 = time.time()
    out_parquet = Path(out_parquet)
    out_parquet.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_parquet.parent / "duckdb_tmp"
    tmp.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    con.execute(f"PRAGMA memory_limit='{memory_limit}'")
    con.execute(f"PRAGMA threads={threads}")
    con.execute("PRAGMA preserve_insertion_order=false")
    con.execute(f"PRAGMA temp_directory='{tmp}'")

    tracks = [(small_vcf, "small")] + ([(sv_vcf, "sv")] if sv_vcf else [])
    union = "\n        UNION ALL\n".join(_track_sql(v, t) for v, t in tracks)

    con.execute(f"""
        COPY (
            WITH ann AS ({union})
            SELECT a.variant_id, a.chrom, a."start", a.ref, a.alt, a.variant_type,
                   a.gene, m.symbol, a.transcript, m.gene_type, m.biotype, m.canonical, m.tsl,
                   a.consequence, a.impact, a.distance, a.exon, a.intron,
                   a.hgvsc, a.hgvsp, a.ensp,
                   a.cdna_position, a.cds_position, a.protein_position,
                   a.amino_acids, a.codons, a.strand
            FROM ann a
            LEFT JOIN read_parquet('{meta}') m ON a.transcript = m.transcript_id
        ) TO '{out_parquet}' (FORMAT parquet)
    """)
    # cheap count from the parquet footer (no COUNT DISTINCT sweep)
    n = con.execute(f"SELECT count(*) FROM read_parquet('{out_parquet}')").fetchone()[0]
    print(f"wrote {out_parquet} ({n:,} variant×transcript rows, "
          f"tracks={[t for _, t in tracks]}) in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    smk = snakemake  # noqa: F821
    sv = getattr(smk.input, "sv", None)
    parse(smk.input.small, sv, smk.input.metadata, smk.output.parquet,
          smk.params.memory_limit, int(smk.params.threads))
