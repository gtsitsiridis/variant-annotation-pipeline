# CLAUDE.md — variant-annotation-pipeline-v2

Living notes for Claude. **Keep this updated** as the pipeline evolves. Last updated: 2026-06-26.

## What this is
The **dataset-agnostic annotation engine**: a generic unique-variant VCF in → annotation
parquets out. NOT GTEx-specific. Consumed across a repo boundary (gtex-benchmark is the dataset
adapter) either **standalone** or as a Snakemake **`module`**. GitLab remote
`gagneurlab/variant-annotation-pipeline-v2`. (⚠️ This **v2** repo superseded an older
`~/projects/variant_annotation`, which was deleted 2026-06-26 — v2 is canonical.)

## Architecture — per-tool, distance-keyed (REDESIGNED 2026-07-16; the "deep tier" is gone)
Each tool has its own `enabled` flag (`fastvep` always-on). ONE top-level **`distance`** drives
fastVEP `--distance`, the deep-tool variant selection, and VEP `--distance`. Uniform output:
**`OUT/distance_<distance>/<tool>.parquet/variant_type={SNV,indel,SV}/`** — each hive-partitioned by
`variant_type`, keyed on **`variant_id`** (`chrom_start_ref_alt` small / `sv_id` SV). `transcript_metadata.parquet`
is distance-independent (top level). Dataset-agnostic: NO freq/svtype — the consumer rejoins those.
- **`fastvep.smk`** (ALWAYS on; only gff3+fasta): fastVEP over the whole set per `{track}` — `small`
  (`input_vcf`) + optional `sv` (`sv_vcf`, gated by `fastvep.include_sv`; **SV is fastVEP-only**).
  `parse_fastvep.py` splits `small`→`variant_type={SNV,indel}` (ref/alt) and `sv`→`SV` →
  `distance_<d>/fastvep.parquet/variant_type=.../`. The base annotation + the `distance` selection source.
- **`vep.smk`** (`vep.enabled`): `select.smk` picks the small set within `distance` of a gene (from
  fastVEP's `distance` annotation) → `analysis_set.vcf.gz`; `chunk_vcf` → chunked VEP (`--gtf`,
  `--distance {config.distance}`) + plugins (LOFTEE/CADD/SpliceAI/PrimateAI/AlphaMissense); `parse_vep.py`
  carries `variant_id`+`variant_type`; `vep_parquet` concats → `distance_<d>/vep.parquet/variant_type=.../`.
  Canonical is NOT merged (join `transcript_metadata` downstream) — so `merge.smk`/`merge_annotations.py`/
  `build_canonical` are **retired/dead**.
- **`e2g.smk`** (`e2g.enabled`): standalone (reads `input_vcf` + the ENCODE-rE2G tables, no selection);
  `annotate_enhancers.py` overlaps variants with enhancers, derives `variant_type` from `variant_id`
  → `distance_<d>/e2g.parquet/variant_type={SNV,indel}/` (small-only, grain variant×target_gene).
- **NMD / AbSplice** — still on the old machinery; **PENDING migration** to `<tool>.parquet/
  variant_type`. The Snakefile keeps them in a `_PENDING` guard (enabling one raises a clear error).

`Snakefile`: `include fastvep.smk` always; `if VEP: include select + vep`; `if E2G: include e2g`;
nmd/absplice pending. `rule all` = fastVEP (+ vep + e2g) partition dirs. `include_sv` is rejected on
any non-fastvep tool.

## Gotchas (learned the hard way)
- **`script:` files must NOT have a line-1 `from __future__ import annotations`** — Snakemake
  prepends a preamble, pushing the future import off line 1 → SyntaxError at runtime. (Bit the
  Tier-0/E2G run 2026-06-26; fixed. Convention: py3.12 conda envs support modern hints natively.)
  build_canonical/parse_vep/merge_annotations/run_absplice all carry the same NB.
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
- **VEP `--gtf`** (no cache): SIFT/PolyPhen/Condel + cached gnomAD AF unavailable; `--canonical`
  inert → reconstruct from GTF (`build_canonical.py`). **NMD-Scanner**: use CSV output (parquet
  writer buggy). Chunked VEP scatter (deeprvat-style) avoids a big-chromosome straggler.
- **Minor debt:** `build_transcript_metadata.py` (Tier 0) overlaps `build_canonical.py` (deep) —
  both derive canonical from the GTF. Kept separate to not disturb the validated deep merge;
  consolidate later.

## Conventions
- Snakemake; per-rule `conda:` envs (`envs/*.yaml`). Run with `--use-conda` (+ `profiles/slurm`
  for the cluster, or `--cores N` local). Tier 0 + DuckDB parse run cheaply; deep VEP is the heavy
  Slurm scatter.
- Config-driven paths; `output_dir`/`input_vcf` absolute. `config.yaml` = template;
  `config.local.yaml` = real paths (gitignored); `config.smoke.yaml` = chr21 check.
- DuckDB for the big explodes/overlaps (parse_fastvep, annotate_enhancers); polars elsewhere.

## Build status (2026-06-26)
- Deep tier (VEP+plugins+NMD+AbSplice, chunked) run end-to-end on the full 30M GTEx set.
- Tier 0 fastVEP + selector funnel + E2G added; tier-gating Snakefile; `--use-conda` envs incl.
  fastvep cargo. Dry-runs pass (Tier 0 / deep / e2g); first full Tier-0+E2G run validating now.
