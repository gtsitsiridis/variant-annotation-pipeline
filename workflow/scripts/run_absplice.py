"""Per-(variant, gene) AbSplice2 scores from the precomputed gagneurlab/absplice result.

We do NOT run AbSplice (its own workflow needs per-tissue SpliceMaps + SpliceAI + GPU); we consume
the precomputed result deeprvat uses: one parquet per gene (~19k `ENSG…_max_preds.parquet`), keyed
on (chrom, start[0-based], ref, alt, gene_id[unversioned]) — SNV-only — carrying `AbSplice_DNA_max`
and `AbSplice2_max` (max across the 49 GTEx tissues).

We inner-join it onto fastVEP's filtered VCF (fastvep/small.vcf.gz) on (chrom, start, ref, alt), so
only the analysis-set variants are scored, and version-map AbSplice's unversioned gene_id to the
gencode versioned gene_id (via the GTF) so combine can join on (variant_id, gene). Output:
variant_id, gene, AbSplice_DNA_max, AbSplice2_max. Indels (absent from the SNV-only source) never
match -> no row (combine leaves them null).

Invoked by Snakemake (snakemake.input.vcf; snakemake.params.result/.gtf/.memory_limit/.threads;
snakemake.output.parquet).
"""

import time
from pathlib import Path

import duckdb


def main(vcf: str, result: str, gtf: str, out_parquet: str,
         memory_limit: str = "16GB", threads: int = 4) -> None:
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

    absp_glob = str(Path(result) / "**" / "*.parquet")
    read_vcf = (f"read_csv('{vcf}', delim='\t', header=false, comment='#', "
                f"all_varchar=true, maximum_line_size=10000000)")
    read_gtf = (f"read_csv('{gtf}', delim='\t', header=false, comment='#', "
                f"all_varchar=true, maximum_line_size=10000000)")

    con.execute(f"""
        COPY (
            WITH v AS (   -- filtered VCF variants: 0-based start = POS-1
                SELECT column2 AS variant_id, column0 AS chrom,
                       CAST(column1 AS BIGINT) - 1 AS start, column3 AS ref, column4 AS alt
                FROM {read_vcf}
            ),
            a AS (        -- precomputed AbSplice2 (max over the 49 tissues), gene_id unversioned
                SELECT CAST(chrom AS VARCHAR) AS chrom, CAST(start AS BIGINT) AS start,
                       CAST(ref AS VARCHAR) AS ref, CAST(alt AS VARCHAR) AS alt,
                       CAST(gene_id AS VARCHAR) AS gene_unversioned,
                       AbSplice_DNA_max, AbSplice2_max
                FROM read_parquet('{absp_glob}', union_by_name=true)
            ),
            gmap AS (     -- gencode unversioned ENSG -> versioned gene_id
                SELECT DISTINCT
                       regexp_extract(column8, 'gene_id "([^"]+)"', 1) AS gene,
                       regexp_replace(regexp_extract(column8, 'gene_id "([^"]+)"', 1), '\\..*$', '')
                           AS gene_unversioned
                FROM {read_gtf} WHERE column2 = 'gene'
            )
            SELECT v.variant_id,
                   coalesce(g.gene, a.gene_unversioned) AS gene,   -- versioned if mapped
                   a.AbSplice_DNA_max, a.AbSplice2_max
            FROM v JOIN a
              ON v.chrom = a.chrom AND v.start = a.start AND v.ref = a.ref AND v.alt = a.alt
            LEFT JOIN gmap g ON a.gene_unversioned = g.gene_unversioned
        ) TO '{out_parquet}' (FORMAT parquet)
    """)
    s = con.execute(f"SELECT count(*), count(DISTINCT variant_id) "
                    f"FROM read_parquet('{out_parquet}')").fetchone()
    print(f"wrote {out_parquet} ({s[0]:,} (variant,gene) rows, {s[1]:,} variants) "
          f"in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    smk = snakemake  # noqa: F821
    p = smk.params
    main(smk.input.vcf, p.result, p.gtf, smk.output.parquet,
         p.memory_limit, int(p.threads))
