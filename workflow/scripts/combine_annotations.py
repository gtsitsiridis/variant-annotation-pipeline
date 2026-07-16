"""Combine the per-method annotation parquets into ONE hive-partitioned annotations.parquet.

Base = the fastVEP transcript table (parts/fastvep.parquet, all variant_types). Optional filters
drop non-matching transcript rows (protein_coding_only -> gene_type == protein_coding;
canonical_only -> canonical). Each enabled `additional` tool is LEFT JOINed on:
  - vep      (variant_id, transcript)  transcript-level  -> vep's plugin/LoF/SpliceAI/Consequence cols
  - e2g      (variant_id, gene)         gene-level        -> e2g_* cols (broadcast across transcripts)
  - absplice (variant_id, gene)         gene-level        -> AbSplice_DNA_max / AbSplice2_max
SV rows (fastVEP-only) never match the additional parts -> those columns stay null. The result is
written PARTITION_BY (variant_type, gene_type, canonical).

Columns from each tool that are join keys or duplicate the fastVEP base are dropped (below); the
tool-specific columns are carried through (e2g's are prefixed to avoid clashes). Which tool parquets
are passed depends on config (combine handles any subset, including none).

Invoked by Snakemake (snakemake.input.fastvep [+ .vep/.e2g/.absplice]; snakemake.params.out_dir/
.protein_coding_only/.canonical_only/.memory_limit/.threads).
"""

import shutil
import time
from pathlib import Path

import duckdb

# tool columns that are join keys or duplicate the fastVEP base -> not carried into the combined table
_VEP_DROP = {"variant_id", "chrom", "start", "ref", "alt", "Feature", "variant_type",
             "Gene", "BIOTYPE", "CANONICAL", "IMPACT", "SYMBOL", "Consequence"}
_E2G_DROP = {"variant_id", "target_gene_id"}
_E2G_PREFIX = "e2g_"
_ABSPLICE_DROP = {"variant_id", "gene"}


def _cols(con, path):
    return [r[0] for r in con.execute(
        f"DESCRIBE SELECT * FROM read_parquet('{path}')").fetchall()]


def main(fastvep: str, out_dir: str, vep=None, e2g=None, absplice=None, *,
         protein_coding_only=False, canonical_only=False,
         memory_limit="64GB", threads=6) -> None:
    t0 = time.time()
    out_dir = Path(out_dir)
    shutil.rmtree(out_dir, ignore_errors=True)          # clean rebuild (drop stale partitions)
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp = out_dir.parent / "duckdb_tmp"
    tmp.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    con.execute(f"PRAGMA memory_limit='{memory_limit}'")
    con.execute(f"PRAGMA threads={threads}")
    con.execute("PRAGMA preserve_insertion_order=false")
    con.execute(f"PRAGMA temp_directory='{tmp}'")

    selects, joins = ["b.*"], []
    if vep:
        extra = [c for c in _cols(con, vep) if c not in _VEP_DROP]
        selects += [f'vep."{c}"' for c in extra]
        joins.append(f"LEFT JOIN read_parquet('{vep}') vep "
                     f'ON b.variant_id = vep.variant_id AND b.transcript = vep."Feature"')
    if e2g:
        extra = [c for c in _cols(con, e2g) if c not in _E2G_DROP]
        selects += [f'e2g."{c}" AS "{_E2G_PREFIX}{c}"' for c in extra]
        joins.append(f"LEFT JOIN read_parquet('{e2g}') e2g "
                     f"ON b.variant_id = e2g.variant_id AND b.gene = e2g.target_gene_id")
    if absplice:
        extra = [c for c in _cols(con, absplice) if c not in _ABSPLICE_DROP]
        selects += [f'absplice."{c}"' for c in extra]
        joins.append(f"LEFT JOIN read_parquet('{absplice}') absplice "
                     f"ON b.variant_id = absplice.variant_id AND b.gene = absplice.gene")

    where = ["TRUE"]
    if protein_coding_only:
        where.append("b.gene_type = 'protein_coding'")
    if canonical_only:
        where.append("b.canonical")

    con.execute(f"""
        COPY (
            SELECT {', '.join(selects)}
            FROM read_parquet('{fastvep}') b
            {' '.join(joins)}
            WHERE {' AND '.join(where)}
        ) TO '{out_dir}'
          (FORMAT parquet, PARTITION_BY (variant_type, gene_type, canonical), OVERWRITE_OR_IGNORE)
    """)
    n = con.execute(
        f"SELECT count(*) FROM read_parquet('{out_dir}/**/*.parquet', hive_partitioning=true)"
    ).fetchone()[0]
    tools = [t for t, p in (("vep", vep), ("e2g", e2g), ("absplice", absplice)) if p]
    print(f"wrote {out_dir} ({n:,} rows; joined {tools or 'none'}; "
          f"protein_coding_only={protein_coding_only} canonical_only={canonical_only}) "
          f"in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    smk = snakemake  # noqa: F821
    inp, p = smk.input, smk.params
    main(inp.fastvep, p.out_dir,
         vep=getattr(inp, "vep", None), e2g=getattr(inp, "e2g", None),
         absplice=getattr(inp, "absplice", None),
         protein_coding_only=bool(p.protein_coding_only),
         canonical_only=bool(p.canonical_only),
         memory_limit=p.memory_limit, threads=int(p.threads))
