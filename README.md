# variant-annotation-pipeline

**Dataset-agnostic** functional annotation of a variant set. Takes a generic
**unique-variant VCF** (`ID = chrom_start_ref_alt`, from any upstream caller — e.g. the
gtex-benchmark `variants` module) and annotates it. Independent of any dataset: it annotates
*variants*, not GTEx, and is consumed across a repo boundary (e.g. by gtex-benchmark's
`evaluation` module) either **standalone** or **as a Snakemake `module`**.

## Annotation model — fastVEP funnel → one `annotations.parquet`

fastVEP is the always-on **driver**. It annotates the whole set (the **base** table), then applies the
pipeline filters (`distance` + `protein_coding_only` + `canonical_only`) to emit a **filtered VCF**
that every enabled `additional` tool re-annotates — so the heavy tools run only on the near-gene /
coding / canonical subset. All per-method outputs are LEFT-joined into ONE `annotations.parquet`,
hive-partitioned by `variant_type` / `gene_type` / `canonical`:

```
output_dir/
  transcript_metadata.parquet                     # gencode GTF -> gene_type/canonical/biotype/… (shared)
  fastvep/
    small.csq.vcf.gz  sv.csq.vcf.gz                # raw fastVEP CSQ (temp; the parse reads them)
    small.vcf.gz (+.tbi)                           # THE fastVEP output VCF: filtered small subset
  parts/  fastvep.parquet vep.parquet e2g.parquet absplice.parquet   # temp per-method
  annotations.parquet/variant_type=<vt>/gene_type=<gt>/canonical=<bool>/   # FINAL (hive)
```

- **fastVEP** (always on; only gencode **GFF3 + FASTA**) — per-transcript consequence/distance over the
  whole set, the base table, AND the funnel: the filtered `fastvep/small.vcf.gz`. The **only SV-capable
  tool** (`fastvep.include_sv`; `sv_vcf` optional). SVs go into the base/annotations but never to the
  additional tools.
- **VEP + plugins** (`additional.vep.enabled`) — LOFTEE/CADD/SpliceAI/PrimateAI/AlphaMissense on the
  filtered VCF (chunked), `--gtf` gencode namespace, own `additional.vep.distance`. Transcript-level →
  joined on (variant_id, transcript).
- **NMD** (`additional.nmd.enabled`) — NMD-Scanner PTC escape prediction (`nmd_escape` + the 5 escape
  rules) on the filtered VCF. Transcript-level → joined on (variant_id, transcript).
- **E2G** (`additional.e2g.enabled`) — ENCODE-rE2G enhancer → target-gene overlap, optional
  `distance_to_tss` cap. Gene-level → joined on (variant_id, gene).
- **AbSplice** (`additional.absplice.enabled`) — join a *precomputed* AbSplice2 result onto the filtered
  VCF (SNV-only). Gene-level → joined on (variant_id, gene). AbSplice is **not** run here.

**Filters + partitioning:** `protein_coding_only` / `canonical_only` are optional — when on they drop
non-matching transcript rows (so those partitions don't appear); either way the output is partitioned by
`variant_type` / `gene_type` / `canonical`. **SV rows** carry the additional-tool columns as **null**
(SVs are fastVEP-only). **NMD** is not yet migrated to this funnel (enabling `additional.nmd` errors).
Dropping the `distance_<d>/` prefix means one annotation set per output dir — widen `distance` for a
larger cis / distal reach rather than keeping multiple.

VEP gene namespace is pinned to the config `gtf` via **`--gtf`** (no cache): SIFT/PolyPhen/Condel +
cached gnomAD AF are unavailable (missense via AlphaMissense/PrimateAI/CADD), and `--canonical`/gene_type
come from `transcript_metadata` (joined in the base).

## Input contract
- `input_vcf` — bgzipped + tabix'd VCF of **unique** variants, `ID = chrom_start_ref_alt`
  (carried through as `variant_id`, the join key).
- `sv_vcf` — **optional** SV VCF (`ID = sv_id`); annotated by fastVEP only. Omit it entirely for a
  no-SV dataset (no `variant_type=SV` partition).
- Reference FASTA + gencode **GFF3** (fastVEP) and **GTF** (VEP + metadata), matching the consumer's
  gene namespace (gtex-benchmark = gencode v34).

## Configuration
`config.yaml` is the template (copy to `config.local.yaml`, or pass config via the `module` block):
- `fastvep: { distance, protein_coding_only, canonical_only, include_sv, scratch, parse_memory_limit,
  parse_threads }` — the always-on driver + its filters.
- `additional: { variants_per_chunk, vep, nmd, e2g, absplice }`. `variants_per_chunk` is the **shared**
  VEP/NMD `chunk_vcf` scatter size. Each tool is `{ enabled: false, … }`: `vep` carries `distance` (its
  own `--distance`) + `plugin_dir`/`fork`/`plugin_data{…}`; `nmd` a `reassign_exons` flag; `e2g` the
  ENCODE-rE2G tables + `distance_to_tss`; `absplice` a precomputed `result`.
- No `include_sv` outside `fastvep` (the Snakefile errors if set on an additional tool).

## Usage — standalone
Each rule declares its own conda env, so run with `--use-conda` (conda/mamba required). fastVEP's
`envs/fastvep.yaml` builds the fastVEP Rust CLI via a pinned cargo post-deploy. Deep-tier reference
data (VEP plugins/LOFTEE/CADD/…) and the E2G tables are installed separately and pointed at in config.

```bash
cp config.yaml config.local.yaml               # edit reference + input paths, flip tool flags
uv run snakemake --configfile config.local.yaml -n                       # dry-run the DAG
uv run snakemake --configfile config.local.yaml --use-conda --cores 8    # fastVEP only (default)
uv run snakemake --configfile config.local.yaml --use-conda --profile profiles/slurm   # with tools enabled
uv run snakemake --configfile config.local.yaml --delete-all-output      # scoped cleanup
```

## Usage — as a Snakemake `module`
A consumer loads this pipeline as a module and drives it from its own Snakefile. **Gotcha: the
`module` directive's `config:` REPLACES the module's config (it does not merge)** — so the parent
assembles the *complete* config. The clean pattern (as gtex-benchmark does):

```python
_VA = config["variant_annotation"]                       # block in the parent's config
_VA_CFG = {
    **{k: v for k, v in _VA.items() if k not in ("repo", "ref")},   # fastvep, additional, …
    "fasta": config["fasta"], "gtf": config["gtf"], "gff3": config["gff3"],
    "input_vcf": f"{VAR}/variants.vcf.gz",
    "output_dir": f"{OUTPUT}/variant_annotation",
}
if HAS_SV:                                                # sv_vcf is OPTIONAL — omit for a no-SV dataset
    _VA_CFG["sv_vcf"] = f"{VAR}/sv_variants.vcf.gz"

module variant_annotation:
    snakefile: github(_VA["repo"], path="Snakefile", branch=_VA["ref"])
    config:    _VA_CFG

use rule * from variant_annotation as va_*               # rules become va_fastvep, va_parse_fastvep, va_combine, …
```

Then depend on the outputs (e.g. `…/variant_annotation/annotations.parquet/variant_type=SNV`)
from the consumer's own rules — the cross-repo DAG wires automatically. Run the parent with
`--use-conda`. Notes: keep the same Snakemake version on both sides; the module's relative
`include:`/`conda:` paths resolve to *this* repo; loading over `github(...)` uses the pushed
commit, so **commit + push module changes before the consumer picks them up**.

## Reference data / envs
VEP + plugins + their data (CADD ~300 GB, SpliceAI, PrimateAI, AlphaMissense, LOFTEE), NMD-Scanner
(conda; use its **CSV** output — the parquet writer is buggy), AbSplice2 (precomputed result dir),
and the ENCODE-rE2G hg38 tables are installed/downloaded separately. fastVEP is built reproducibly
inside `envs/fastvep.yaml` (rust ≥1.88 + cargo post-deploy, pinned).
