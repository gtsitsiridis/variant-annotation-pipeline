# CLAUDE.md — variant-annotation-pipeline-v2

Living notes for Claude. **Keep this updated** as the pipeline evolves. Last updated: 2026-07-16.

## What this is
The **dataset-agnostic annotation engine**: a generic unique-variant VCF in → annotation
parquets out. NOT GTEx-specific. Consumed across a repo boundary (gtex-benchmark is the dataset
adapter) either **standalone** or as a Snakemake **`module`**. GitLab remote
`gagneurlab/variant-annotation-pipeline-v2`. (⚠️ This **v2** repo superseded an older
`~/projects/variant_annotation`, which was deleted 2026-06-26 — v2 is canonical.)

## Architecture — fastVEP funnel + combined annotations.parquet (REWORKED 2026-07-16)
fastVEP is the always-on **driver/funnel**. It annotates the whole set, becomes the **base** table,
AND emits a **filtered VCF** that every `additional` tool re-annotates; combine LEFT-joins them into
ONE `annotations.parquet`. Config = `fastvep:` (the driver + its filters) + `additional:` (vep/e2g/
absplice, each `enabled`). NO `distance_<d>/` prefix. Output layout:
```
OUT/ transcript_metadata.parquet          # gencode GTF -> gene_type/canonical/biotype/tsl/symbol (shared)
     fastvep/ small.csq.vcf.gz sv.csq.vcf.gz (temp)   small.vcf.gz (persisted filtered VCF)
     parts/  fastvep.parquet vep.parquet e2g.parquet absplice.parquet   (temp per-method)
     annotations.parquet/variant_type=<vt>/gene_type=<gt>/canonical=<bool>/   (FINAL, hive)
```
- **`fastvep.smk`** (ALWAYS): `rule fastvep` annotates each `{track}` (`small`=input_vcf, `sv`=sv_vcf if
  `fastvep.include_sv`) → temp **bgzipped** `{track}.csq.vcf.gz` (`fastvep -o /dev/stdout | bgzip`).
  `parse_fastvep.py` UNIONs the tracks → ONE flat `parts/fastvep.parquet` (variant×transcript, a
  `variant_type` col, joined to transcript_metadata for gene_type/canonical). `select_analysis_set.py`
  filters the small rows (distance + protein_coding_only + canonical_only) → ids → `bcftools view` →
  **`fastvep/small.vcf.gz`** (the funnel output; only built when an additional tool needs it).
- **`vep.smk`** (`additional.vep.enabled`): `chunk_vcf` splits `fastvep/small.vcf.gz` → chunked VEP
  (`--gtf`, `--distance {additional.vep.distance}`) + plugins (LOFTEE/CADD/SpliceAI/PrimateAI/AlphaMissense);
  `parse_vep.py` (variant_id + `Feature`=transcript + plugin/LoF/SpliceAI/Consequence cols) → concat →
  `parts/vep.parquet`.
- **`nmd.smk`** (`additional.nmd.enabled`): `nmd-scanner` on `fastvep/small.vcf.gz` → CSV; `parse_nmd`
  rebuilds `variant_id` (chrom_start_ref_alt) + transcript + the escape flags → `parts/nmd.parquet`
  (transcript-level). Single job on the filtered subset (chunk like VEP if it straggles).
- **`e2g.smk`** (`additional.e2g.enabled`): `annotate_enhancers.py` overlaps `fastvep/small.vcf.gz` with
  ENCODE-rE2G enhancers, optional `distance_to_tss` cap → `parts/e2g.parquet` (variant×target_gene, gene-level).
- **`absplice.smk`** (`additional.absplice.enabled`): `run_absplice.py` joins the precomputed AbSplice2
  result onto `fastvep/small.vcf.gz` on (chrom,start,ref,alt), version-maps the gene → `parts/absplice.parquet`
  (variant×gene, gene-level, SNV-only).
- **`combine.smk`** (ALWAYS): `combine_annotations.py` (DuckDB) = base `parts/fastvep.parquet` (+ filters)
  LEFT JOIN vep on (variant_id, transcript) + e2g/absplice on (variant_id, gene) → `annotations.parquet`
  PARTITION_BY (variant_type, gene_type, canonical). SV rows: additional cols null (SVs are fastVEP-only).

`Snakefile`: `include fastvep.smk` always; `if additional.<tool>.enabled: include vep/nmd/e2g/absplice`;
`include combine.smk` always. `rule all` = the `annotations.parquet/variant_type=…` dirs. `include_sv`
rejected on any additional tool. Retired: `select.smk`, `merge.smk`/`merge_annotations.py`,
`build_canonical.py`, the `distance_<d>/<tool>.parquet` per-tool layout.

## Gotchas (learned the hard way)
- **`script:` files must NOT have a line-1 `from __future__ import annotations`** — Snakemake
  prepends a preamble, pushing the future import off line 1 → SyntaxError at runtime. (Bit the
  run 2026-06-26; fixed. Convention: py3.12 conda envs support modern hints natively.)
  parse_fastvep/parse_vep/annotate_enhancers/run_absplice/combine_annotations all carry the same NB.
- **Module `config:` REPLACES** the module's config (does not merge with `configfile:`). The
  parent must pass the COMPLETE config (see README "as a module"). This repo's own `configfile:`
  is inert under module loading.
- **fastVEP = conda+cargo** (`envs/fastvep.yaml` + `fastvep.post-deploy.sh`): rust ≥1.88 builds
  the CLI (Huang-lab/fastVEP, pinned rev). `--use-conda` triggers it. The annotate step uses a
  `.fastvep.cache` next to the gff3 (built cold on first use).
- **E2G data** (`/s/raw/e2g/`): a star schema — `enhancer_gene_predictions/chr*.parquet`
  (enhancer_id→target_gene + rE2G score) + `enhancers.parquet` (id→chrom/start/end/class). hg38;
  `chromosome` is **UNprefixed** ("21"); `target_gene_id` is **unversioned** → mapped to gencode
  versioned via the gtf. Models: ENCODE-rE2G (default) + scE2G. Default keeps `genic`+`intergenic`
  classes (drops promoter), score≥0.5, all cell_types (heavy genome-wide — restrict for speed).
- **VEP `--gtf`** (no cache): SIFT/PolyPhen/Condel + cached gnomAD AF unavailable; `--canonical` inert
  → canonical/gene_type come from `transcript_metadata` (joined in `parse_fastvep`, the combine base).
  Chunked VEP scatter (deeprvat-style) avoids a big-chromosome straggler.
- **combine join keys:** vep on `(variant_id, transcript==vep.Feature)`; e2g/absplice on `(variant_id,
  gene)` (gene-level → broadcast across the gene's transcript rows). `combine_annotations._*_DROP` lists
  the tool columns dropped as join-keys / base-duplicates; e2g cols get an `e2g_` prefix.
- **partition boolean:** `annotations.parquet` partitions on `canonical` (bool) → `canonical=true/false`
  dirs; hive readers may surface it as a string — cast/compare as needed.

## Conventions
- Snakemake; per-rule `conda:` envs (`envs/*.yaml`). Run with `--use-conda` (+ `profiles/slurm`
  for the cluster, or `--cores N` local). Tier 0 + DuckDB parse run cheaply; deep VEP is the heavy
  Slurm scatter.
- Config-driven paths; `output_dir`/`input_vcf` absolute. `config.yaml` = template;
  `config.local.yaml` = real paths (gitignored); `config.smoke.yaml` = chr21 check.
- DuckDB for the big explodes/overlaps (parse_fastvep, annotate_enhancers); polars elsewhere.

## Build status (2026-07-16)
- Reworked to the fastVEP-funnel + combined `annotations.parquet` model (this session). fastvep base →
  filtered VCF → vep/e2g/absplice → combine (PARTITION_BY variant_type/gene_type/canonical).
- Dry-runs pass: fastvep-only (5 jobs) and all-tools (11 jobs). `combine_annotations.py` +
  `parse_fastvep.py` synthetic-tested (joins, SV-null additional cols, hive round-trip, filters).
- All four additional tools (vep/nmd/e2g/absplice) migrated to the funnel (transcript-level vep/nmd →
  join on (variant_id, transcript); gene-level e2g/absplice → join on (variant_id, gene)).
- NOT yet run on the cluster with the real conda tools (VEP/NMD/AbSplice reference data + envs).
