# variant-annotation-pipeline

**Dataset-agnostic** functional annotation of a variant set. Takes a generic
**unique-variant VCF** (`ID = chrom_start_ref_alt`, from any upstream caller — e.g. the
gtex-benchmark `variants` module) and annotates it. Independent of any dataset: it annotates
*variants*, not GTEx, and is consumed across a repo boundary (e.g. by gtex-benchmark's
`evaluation` module) either **standalone** or **as a Snakemake `module`**.

## Tiered funnel

```
input_vcf (+ optional sv_vcf)
   │
   ├─ Tier 0  fastVEP        broad per-transcript consequence/distance over the WHOLE set,
   │                         per track (small + SV) — cheap, ALWAYS on, only gff3+fasta
   │                            → <track>/basic_annotations.parquet  (+ transcript_metadata.parquet)
   │            │ select.smk (deep.window) from the small track
   │  deep   VEP + plugins   LOFTEE/CADD/SpliceAI/PrimateAI/AlphaMissense (+ NMD, AbSplice),
   │                         chunked, on the narrowed analysis_set    → annotations.parquet
   │
   └─ E2G    ENCODE-rE2G     per-variant enhancer → target gene + rE2G score (hg38)
                                                                      → e2g/e2g.parquet
```

- **Tier 0 (fastVEP)** — always on; needs only the gencode **GFF3 + FASTA** (no plugin
  reference data). Runs per `{track}`: `small` (`input_vcf`, ID=chrom_start_ref_alt) and, if
  `sv_vcf` is set, `sv` (ID=sv_id; fastVEP is ID-agnostic). It's both the broad annotation
  layer and the basis for selecting the deep-tier window.
- **deep tier** (`deep.enabled`) — `select.smk` narrows the small track to variants within
  `deep.window` bp of a transcript → `analysis_set.vcf.gz`; the VEP scatter (`chunk_vcf`)
  reads *that*, so the heavy plugins only run where they can score. VEP `--gtf` (gencode
  namespace) + CADD/SpliceAI/PrimateAI/AlphaMissense/LOFTEE, NMD-Scanner, and a precomputed
  AbSplice2 lookup are merged into `annotations.parquet`. `nmd`/`absplice` are independently
  gated.
- **E2G tier** (`e2g.enabled`) — overlaps the variants with **ENCODE-rE2G** hg38 enhancer
  elements (the on-disk star-schema export: per-chrom predictions + `enhancers.parquet`),
  giving each variant its predicted target gene + rE2G score + distance to TSS. Kept a
  **standalone** parquet (grain = variant × target_gene; the consumer joins it on
  `variant_id`) — *not* folded into the transcript-level `annotations.parquet`.

VEP gene namespace is pinned to the config `gtf` via **`--gtf`** (no cache), so
`Gene`/`Feature` match the consumer (gencode v34). Trade-off: SIFT/PolyPhen/Condel + cached
gnomAD AF are unavailable (cache-only); missense via AlphaMissense/PrimateAI/CADD. VEP's
`--canonical` is inert in `--gtf` mode, so `CANONICAL` is reconstructed from the GTF
(`build_canonical.py`: MANE_Select else best `appris_principal_N`).

## Input contract
- `input_vcf` — bgzipped + tabix'd VCF of **unique** variants, `ID = chrom_start_ref_alt`
  (carried through as VEP `Uploaded_variation`; the join key).
- `sv_vcf` (optional) — SV VCF for Tier 0 only (`ID = sv_id`).
- Reference FASTA + gencode **GFF3** (Tier 0) and **GTF** (VEP/canonical/metadata), matching
  the consumer's gene namespace (gtex-benchmark = gencode v34).

## Outputs (under `output_dir`)
| File | Tier | Grain |
|---|---|---|
| `<track>/basic_annotations.parquet` | 0 | variant × transcript (keyed on `variant_id`) |
| `transcript_metadata.parquet` | 0 | transcript (symbol/biotype/canonical/tsl) |
| `annotations.parquet` | deep | variant × transcript (VEP + plugins [+ NMD + AbSplice]) |
| `e2g/e2g.parquet` | E2G | variant × target_gene (in-enhancer + rE2G score + dist) |

## Configuration
`config.yaml` is the template (copy to `config.local.yaml` and edit, or pass config via the
`module` block — see below). Tier flags: `deep.enabled` / `nmd.enabled` / `absplice.enabled`
/ `e2g.enabled` (Tier 0 always runs). Key blocks: `fastvep` (distance, scratch), `deep`
(window), `vep` + `plugin_data` (deep reference data), `e2g` (predictions/enhancers dirs,
model, score_threshold, classes, cell_types).

## Usage — standalone
Each rule declares its own conda env, so run with `--use-conda` (conda/mamba required).
Tier 0's `envs/fastvep.yaml` builds the fastVEP Rust CLI via a cargo post-deploy (pinned).
Deep-tier reference data (VEP plugins/LOFTEE/CADD/…) and the E2G tables are installed
separately and pointed at in the config.

```bash
# configure paths
cp config.yaml config.local.yaml      # edit reference + input paths, flip tier flags

# dry-run the DAG
uv run snakemake --configfile config.local.yaml -n

# Tier 0 only (default; just gff3+fasta) — local
uv run snakemake --configfile config.local.yaml --use-conda --cores 8

# full run (deep + E2G enabled in the config) on Slurm
uv run snakemake --configfile config.local.yaml --use-conda --profile profiles/slurm

# scoped cleanup (no rm -rf)
uv run snakemake --configfile config.local.yaml --delete-all-output
```

## Usage — as a Snakemake `module`
A consumer loads this pipeline as a module and drives it from its own Snakefile. **Gotcha:
the `module` directive's `config:` REPLACES the module's config (it does not merge with the
module's own `configfile:`)** — so the parent must assemble this pipeline's *complete* config.
The clean pattern (as gtex-benchmark does): keep a `variant_annotation:` block in the parent
config holding the pipeline-specific keys, then add the shared refs + seam paths:

```python
_VA = config["variant_annotation"]                       # block in the parent's config
_VA_CFG = {
    **{k: v for k, v in _VA.items() if k != "repo"},     # chromosomes, fastvep, deep, vep,
                                                          # plugin_data, nmd, absplice, e2g, …
    "fasta": config["fasta"], "gtf": config["gtf"], "gff3": config["gff3"],
    "input_vcf": f"{VAR}/variants.vcf.gz",
    "sv_vcf":    f"{VAR}/sv_variants.vcf.gz",
    "output_dir": f"{OUTPUT}/variant_annotation",
}

module variant_annotation:
    snakefile: f"{_VA['repo']}/Snakefile"
    config:    _VA_CFG

use rule * from variant_annotation as va_*               # rules become va_fastvep, va_vep, …
```

Then depend on the outputs (e.g. `…/variant_annotation/small/basic_annotations.parquet`) from
the consumer's own rules — the cross-repo DAG wires automatically (the consumer's VCF →
`va_fastvep` → …). Run the parent with `--use-conda` so the module's per-rule conda envs
build. Notes: keep the same Snakemake version on both sides; the module's relative
`include:`/`conda:` paths resolve to *this* repo; this repo's own `configfile:` is inert under
module loading.

## Reference data / envs
VEP + plugins + their data (CADD ~300 GB, SpliceAI, PrimateAI, AlphaMissense, LOFTEE), NMD-
Scanner (conda; use its **CSV** output — the parquet writer is buggy), AbSplice2 (precomputed
result dir), and the ENCODE-rE2G hg38 tables are installed/downloaded separately. fastVEP is
built reproducibly inside `envs/fastvep.yaml` (rust ≥1.88 + cargo post-deploy, pinned).

> Status: deep tier (VEP+plugins+NMD+AbSplice, chunked) **run end-to-end on the full 30M-variant
> GTEx set**. Tier 0 fastVEP + E2G added 2026-06-26 (tiered funnel; opt-in deep/E2G).
