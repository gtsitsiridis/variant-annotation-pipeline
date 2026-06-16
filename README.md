# variant_annotation

**Dataset-agnostic** functional annotation of a variant set. Takes a generic
**unique-variant VCF** (e.g. `variants.vcf.gz` from any upstream caller — in particular
the gtex-benchmark `variants` module) and produces a single per-variant annotation table.

This repo is deliberately independent of any dataset: it annotates *variants*, not GTEx.
It is consumed across a repo boundary (e.g. by gtex-benchmark's `evaluation` module, which
needs LOFTEE pLoF / AbSplice / missense scores).

Modeled on the DeepRVAT annotation pipeline (PMBio/deeprvat, `pipelines/annotations*`).

## Annotations

Built in tiers; each tier is an independently runnable Snakemake target.

### Tier 1 — VEP + plugins (one VEP pass)
| Annotation | Source | Output fields |
|---|---|---|
| Consequence / IMPACT / gene / transcript / biotype / canonical | VEP core | one-hot `Consequence_*`, `IMPACT`, `Gene`, `Feature`, … |
| Missense deleteriousness | VEP `--sift --polyphen`, **Condel**, **PrimateAI**, **AlphaMissense** | `sift_score`, `polyphen_score`, `Condel`, `PrimateAI_score`, `am_pathogenicity` |
| Deleteriousness (genome-wide) | **CADD** plugin | `CADD_raw`, `CADD_PHRED` |
| Splicing | **SpliceAI** plugin | `SpliceAI_delta_score` (max of DS_AG/AL/DG/DL) |
| Loss-of-function confidence | **LOFTEE** plugin | `LoF` (HC/LC), `LoF_filter`, `LoF_flags` |
| Allele frequency | VEP gnomAD (`--af_gnomadg`) | `gnomADg_AF`, derived `MAF` |

### Tier 2 — AbSplice (tissue-specific splicing)
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
```bash
# 1. configure paths (VEP cache + plugins, plugin data, reference, input VCF)
cp config.yaml config.local.yaml   # edit
# 2. dry-run
snakemake -n --configfile config.local.yaml
# 3. run (Slurm; tune profile)
snakemake --configfile config.local.yaml --cores 16
```

VEP, its cache/plugins, and the plugin reference data (CADD ~300 GB, SpliceAI, PrimateAI,
AlphaMissense, LOFTEE data) must be installed/downloaded separately — see `config.yaml`
and `envs/`. AbSplice has its own setup (gagneurlab/absplice).

> Status: scaffold. The Snakemake rules and parsing scripts are wired but **not yet
> smoke-tested end-to-end** (needs a VEP install + reference data).
>>>>>>> 7f965a4 (Scaffold variant_annotation: dataset-agnostic VEP + plugins + AbSplice (Tier 1+2))
