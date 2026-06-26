# CLAUDE.md — variant-annotation-pipeline-v2

Living notes for Claude. **Keep this updated** as the pipeline evolves. Last updated: 2026-06-26.

## What this is
The **dataset-agnostic annotation engine**: a generic unique-variant VCF in → annotation
parquets out. NOT GTEx-specific. Consumed across a repo boundary (gtex-benchmark is the dataset
adapter) either **standalone** or as a Snakemake **`module`**. GitLab remote
`gagneurlab/variant-annotation-pipeline-v2`. (⚠️ This **v2** repo superseded an older
`~/projects/variant_annotation`, which was deleted 2026-06-26 — v2 is canonical.)

## Architecture — tiered funnel (Snakefile gates each tier by a config flag)
- **Tier 0 `fastvep.smk`** (ALWAYS on; only gff3+fasta): fastVEP over the WHOLE set, per
  `{track}` — `small` (`input_vcf`, ID=chrom_start_ref_alt) + optional `sv` (`sv_vcf`, ID=sv_id;
  fastVEP is ID-agnostic). `build_transcript_metadata.py` (symbol/biotype/canonical/tsl) +
  `parse_fastvep.py` (CSQ explode, DuckDB) → `<track>/basic_annotations.parquet` keyed on
  `variant_id` (dataset-agnostic: NO freq/variant_type/svtype/`end` — the consumer rejoins those).
- **deep** (`deep.enabled`): `select.smk` narrows the small track to `deep.window` bp of a
  transcript → `analysis_set.vcf.gz`; the checkpoint `chunk_vcf` reads THAT (not full input), so
  chunked VEP (`--gtf`) + plugins (LOFTEE/CADD/SpliceAI/PrimateAI/AlphaMissense) + NMD + AbSplice
  run only on the narrowed set → `annotations.parquet` (merge spines on chunked vep + canonical).
  `nmd`/`absplice` independently gated.
- **E2G** (`e2g.enabled`): `annotate_enhancers.py` overlaps variants with ENCODE-rE2G hg38
  enhancers → `e2g/e2g.parquet`, grain variant×target_gene. STANDALONE (joined on `variant_id`,
  NOT merged into the transcript-level annotations — different grain).

`Snakefile`: `include fastvep.smk` always; `if DEEP: include select/vep/(nmd)/merge`;
`if E2G: include e2g`. `rule all` = per-track basic_annotations (+ annotations.parquet if deep,
+ e2g/e2g.parquet if e2g).

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
