# variant-annotation-pipeline

**Dataset-agnostic** functional annotation of a variant set. Takes a generic
**unique-variant VCF** (`ID = chrom_start_ref_alt`, from any upstream caller — e.g. the
gtex-benchmark `variants` module) and annotates it. Independent of any dataset: it annotates
*variants*, not GTEx, and is consumed across a repo boundary (e.g. by gtex-benchmark's
`evaluation` module) either **standalone** or **as a Snakemake `module`**.

## Per-tool annotation

fastVEP is the always-on base; every other tool is an independent opt-in flag (the old
"deep tier" umbrella is gone). A **single `distance`** (the cis window) drives fastVEP's
`--distance`, the deep-tool variant *selection* (via fastVEP's `distance` annotation column),
and VEP's `--distance`. Output is uniform and **distance-keyed**:

```
output_dir/
  distance_<distance>/
    fastvep.parquet/variant_type={SNV,indel,SV}/     (ALWAYS — the base annotation)
    vep.parquet/variant_type={SNV,indel}/            (vep.enabled)
    nmd.parquet/variant_type={SNV,indel}/            (nmd.enabled)
    absplice.parquet/variant_type={SNV,indel}/       (absplice.enabled — VEP-dependent, SNV-only)
    e2g.parquet/variant_type={SNV,indel}/            (e2g.enabled)
  transcript_metadata.parquet                        (distance-independent, shared)
```

Every per-tool parquet is hive-partitioned by `variant_type` and keyed on **`variant_id`**
(`chrom_start_ref_alt` for small variants, `sv_id` for SVs) — the single join key. Keying the
output by distance lets different-window runs coexist.

- **fastVEP** (always on; only needs gencode **GFF3 + FASTA**) — broad per-transcript
  consequence/distance over the whole set. The **only SV-capable tool** (`fastvep.include_sv`;
  `sv_vcf` is optional). `small` → `variant_type={SNV,indel}` (split by ref/alt); `sv` →
  `variant_type=SV`. It is both the base annotation and the source of the `distance` used to
  select each deep tool's variant set.
- **VEP + plugins** (`vep.enabled`) — LOFTEE/CADD/SpliceAI/PrimateAI/AlphaMissense on the
  fastVEP-selected small set (chunked), `--gtf` gencode namespace. Small-variant only.
- **NMD** (`nmd.enabled`) — NMD-Scanner escape prediction on that set. Small-variant only.
- **AbSplice** (`absplice.enabled`) — join of a *precomputed* AbSplice2 result onto VEP's
  variant×gene keys (SNV-only; depends on VEP). AbSplice is **not** run here.
- **E2G** (`e2g.enabled`) — ENCODE-rE2G hg38 enhancer → target-gene overlap. Small-variant only.

> **Migration status (2026-07-16):** fastVEP, **VEP** and **e2g** are on the layout above
> (`variant_id`, distance-keyed). **NMD / AbSplice are pending** — gated (enabling one raises an
> error) until migrated. **SV annotation is fastVEP-only** (setting `include_sv` on any other tool
> is rejected).

VEP gene namespace is pinned to the config `gtf` via **`--gtf`** (no cache), so `Gene`/`Feature`
match the consumer (gencode v34). Trade-off: SIFT/PolyPhen/Condel + cached gnomAD AF are
unavailable (cache-only); missense via AlphaMissense/PrimateAI/CADD. VEP's `--canonical` is
inert in `--gtf` mode, so `CANONICAL` is reconstructed from the GTF (MANE_Select else best
`appris_principal_N`).

## Input contract
- `input_vcf` — bgzipped + tabix'd VCF of **unique** variants, `ID = chrom_start_ref_alt`
  (carried through as `variant_id`, the join key).
- `sv_vcf` — **optional** SV VCF (`ID = sv_id`); annotated by fastVEP only. Omit it entirely for
  a no-SV dataset (the pipeline just skips the `variant_type=SV` partition).
- Reference FASTA + gencode **GFF3** (fastVEP) and **GTF** (VEP/canonical/metadata), matching the
  consumer's gene namespace (gtex-benchmark = gencode v34).

## Configuration
`config.yaml` is the template (copy to `config.local.yaml`, or pass config via the `module`
block). Keys:
- `distance` — the one cis window (fastVEP `--distance` + deep selection + VEP `--distance`); also
  the output-dir prefix.
- `fastvep: { include_sv, scratch, parse_memory_limit, parse_threads }` — always on.
- `vep / nmd / absplice / e2g: { enabled: false, … }` — per-tool opt-in. `vep` also carries
  `plugin_dir`/`fork`/`variants_per_chunk`/`plugin_data{…}`; `absplice` a precomputed `result`;
  `e2g` the ENCODE-rE2G tables.
- No per-tool `distance`; no `include_sv` outside `fastvep` (the Snakefile errors if set elsewhere).

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
    **{k: v for k, v in _VA.items() if k not in ("repo", "ref")},   # distance, fastvep, vep, nmd, absplice, e2g, …
    "fasta": config["fasta"], "gtf": config["gtf"], "gff3": config["gff3"],
    "input_vcf": f"{VAR}/variants.vcf.gz",
    "output_dir": f"{OUTPUT}/variant_annotation",
}
if HAS_SV:                                                # sv_vcf is OPTIONAL — omit for a no-SV dataset
    _VA_CFG["sv_vcf"] = f"{VAR}/sv_variants.vcf.gz"

module variant_annotation:
    snakefile: github(_VA["repo"], path="Snakefile", branch=_VA["ref"])
    config:    _VA_CFG

use rule * from variant_annotation as va_*               # rules become va_fastvep, va_parse_fastvep_small, …
```

Then depend on the outputs (e.g. `…/variant_annotation/distance_<d>/fastvep.parquet/variant_type=SNV`)
from the consumer's own rules — the cross-repo DAG wires automatically. Run the parent with
`--use-conda`. Notes: keep the same Snakemake version on both sides; the module's relative
`include:`/`conda:` paths resolve to *this* repo; loading over `github(...)` uses the pushed
commit, so **commit + push module changes before the consumer picks them up**.

## Reference data / envs
VEP + plugins + their data (CADD ~300 GB, SpliceAI, PrimateAI, AlphaMissense, LOFTEE), NMD-Scanner
(conda; use its **CSV** output — the parquet writer is buggy), AbSplice2 (precomputed result dir),
and the ENCODE-rE2G hg38 tables are installed/downloaded separately. fastVEP is built reproducibly
inside `envs/fastvep.yaml` (rust ≥1.88 + cargo post-deploy, pinned).
