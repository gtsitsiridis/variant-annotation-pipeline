# variant_annotation

**Dataset-agnostic** functional annotation of a variant set. Takes a generic
**unique-variant VCF** (e.g. `variants.vcf.gz` from any upstream caller — in particular
the gtex-benchmark `variants` module) and produces a single per-variant annotation table.

This repo is deliberately independent of any dataset: it annotates *variants*, not GTEx.
It is consumed across a repo boundary (e.g. by gtex-benchmark's `evaluation` module, which
needs LOFTEE pLoF / AbSplice / missense scores).

Modeled on the DeepRVAT annotation pipeline (PMBio/deeprvat, `pipelines/annotations*`).

## Annotations

Each annotation group is an independently runnable Snakemake target.

### VEP + plugins (one VEP pass)
| Annotation | Source | Output fields |
|---|---|---|
| Consequence / IMPACT / gene / transcript / biotype / canonical | VEP core (`--gtf` gencode v34) | one-hot `Consequence_*`, `IMPACT`, `Gene`, `Feature`, exon/intron, … |
| Missense deleteriousness | **PrimateAI**, **AlphaMissense** | `PrimateAI_score`, `am_pathogenicity` |
| Deleteriousness (genome-wide) | **CADD** plugin | `CADD_raw`, `CADD_PHRED` |
| Splicing | **SpliceAI** plugin | `SpliceAI_delta_score` (max of DS_AG/AL/DG/DL) |
| Loss-of-function confidence | **LOFTEE** plugin | `LoF` (HC/LC), `LoF_filter`, `LoF_flags` |
| **NMD escape** (PTC fate) | **NMD-Scanner** (gagneurlab) | `nmd_escape` + 5 escape rules (last-exon, 50nt-penultimate, long-exon, start-proximal, single-exon), `alt_is_premature`, `start_loss`, `stop_loss` |

Gene namespace is pinned to the config `gtf` (gencode v34) via VEP **`--gtf`** (custom
annotation), so `Gene`/`Feature` match the consumer. Trade-off of `--gtf` (no cache):
**SIFT/PolyPhen/Condel and cached gnomAD AF are unavailable** (cache-only) — missense
deleteriousness comes from AlphaMissense/PrimateAI/CADD, and allele frequency from the
source variant set (or add a `--custom` gnomAD VCF). Verified live: VEP 113 `--gtf` on a
chr21 subset yields gencode-versioned `Gene`/`Feature` (`ENSG…​.N`/`ENST…​.N`) parsed
straight into the variant key.

### AbSplice (tissue-specific splicing)
`AbSplice_DNA` per tissue (gagneurlab/absplice). Pairs naturally with tissue-specific
downstream analyses (e.g. GTEx tissues).

## Input contract
- `input_vcf` — bgzipped, tabix-indexed VCF of **unique** variants, with `ID =
  chrom_start_ref_alt` (the gtex-benchmark `variants` module writes exactly this). The ID
  is the join key carried through as VEP `Uploaded_variation`.
- Reference FASTA + gencode GTF matching the consumer's gene namespace (gencode v34 for
  gtex-benchmark).

## Output
`output/annotations.parquet` — transcript-level VEP rows (one per variant×transcript)
joined with the variant-level scores, keyed on `chrom,start,ref,alt` (+ `Feature` for the
transcript-level fields). Per-chromosome intermediates under `output/`.

## Run
Runs go through a Snakemake **Slurm profile** (`profiles/slurm/`, mirroring gtex-benchmark),
which submits each job to the `standard` partition. VEP/bcftools/tabix come from the
`vep_v113` conda env put on `PATH` at submit time (inherited by jobs); for a local run
override with `--executor local --cores N`. Dry-runs (`-n`) just resolve the DAG.
```bash
# 1. configure paths (VEP plugins, plugin data, reference, input VCF)
cp config.yaml config.local.yaml   # edit
export PATH=/opt/modules/i12g/anaconda/envs/vep_v113/bin:$PATH
PROFILE="--profile profiles/slurm"

# 2. dry-run (resolve the DAG)
uv run snakemake $PROFILE --configfile config.local.yaml -n

# 3. run the full pipeline on Slurm
uv run snakemake $PROFILE --configfile config.local.yaml

# 4. quick chr21-only check (VEP only; NMD/AbSplice off — see config.smoke.yaml)
uv run snakemake $PROFILE --configfile config.smoke.yaml

# 5. delete only the workflow's declared outputs (scoped; no rm -rf)
uv run snakemake $PROFILE --configfile config.local.yaml --delete-all-output
```

VEP, its cache/plugins, and the plugin reference data (CADD ~300 GB, SpliceAI, PrimateAI,
AlphaMissense, LOFTEE data) must be installed/downloaded separately — see `config.yaml`
and `envs/`. AbSplice has its own setup (gagneurlab/absplice).

> Status: **run end-to-end on the full 30M-variant GTEx set.** `annotations.parquet` =
> 30,422,012 unique variants × transcripts (122M rows, chr1–22+X) with CADD, SpliceAI,
> AlphaMissense, LOFTEE, one-hot consequences, NMD-escape, and AbSplice all populated.
> - **VEP + plugins**: `chunk_vcf` → `vep` (gencode v34 `--gtf` + CADD/SpliceAI/AlphaMissense/
>   LOFTEE, `vep_v113` env) → `parse_vep` → `merge_vep`, on Slurm via `profiles/slurm`.
> - **Chunk-based scatter** (deeprvat-style): the input VCF is split into
>   `vep.variants_per_chunk` chunks (checkpoint `chunk_vcf`), one VEP job each — even
>   parallelism, no big-chromosome straggler. On a cluster, launch the orchestrator itself
>   as a job (`sbatch --wrap '… snakemake --profile profiles/slurm …'`) so it outlives the
>   submitting shell.
> - **NMD-Scanner**: `nmd_escape` + 5 escape-rule flags, joined onto 100% of
>   stop_gained rows (28,787 NMD-escaping PTCs). Needs `nmd-scanner` in a conda env (its
>   `pyranges`/`pysam` deps; use its **CSV** output — the parquet writer is buggy).
> - **AbSplice**: `AbSplice_DNA_max` / `AbSplice2_max` from a precomputed AbSplice2 result,
>   joined per (variant, gene); SNV-only (indels null).
