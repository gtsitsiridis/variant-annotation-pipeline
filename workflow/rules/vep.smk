"""Ensembl VEP + plugins, scattered over fixed-size variant chunks.

Runs on fastVEP's filtered output VCF (`fastvep/small.vcf.gz`) — the near-gene / coding /
canonical subset — so the heavy plugins only score where they can. Gene namespace is pinned to
the config `gtf` (gencode v34) via VEP `--gtf` custom annotation (not a prebuilt cache), keeping
`Gene`/`Feature` in the consumer's namespace. VEP `--gtf` needs a sorted/bgzipped/tabixed GTF
(`prepare_gtf`).

Scatter: the filtered VCF is split into fixed-size chunks of `additional.vep.variants_per_chunk`
variants (checkpoint `chunk_vcf`), one VEP job per chunk (deeprvat-style — no big-chromosome
straggler). Output is per-transcript (no --pick / --per_gene / --most_severe); every transcript
consequence is kept. The per-chunk parquets are concatenated into `parts/vep.parquet`, which
combine joins onto the fastVEP base on (variant_id, transcript).
"""
import os
import sys

VEP_CFG = config["additional"]["vep"]
VEP_PD = VEP_CFG["plugin_data"]
VEP_DISTANCE = VEP_CFG.get("distance", FASTVEP_DISTANCE)   # VEP's own --distance (native)
REF = OUT / "ref"
CHUNKS = OUT / "vcf"                     # chunk_vcf writes chunk_NNNNN.vcf.gz here
CHUNK_SIZE = int(VEP_CFG.get("variants_per_chunk", 500_000))
# Optional chromosome filter applied before chunking ([] / unset = the whole filtered VCF).
_REGION = ("-r " + ",".join(config["chromosomes"])) if config.get("chromosomes") else ""


def chunk_ids():
    """Resolve the dynamic chunk wildcards produced by the chunk_vcf checkpoint."""
    cdir = checkpoints.chunk_vcf.get().output.chunkdir
    return sorted(glob_wildcards(os.path.join(cdir, "chunk_{chunk}.vcf.gz")).chunk)


def _have(*paths) -> bool:
    """True only if every referenced data file is configured and present on disk."""
    return all(p and os.path.exists(p) for p in paths)


def _plugin_args() -> str:
    """Assemble the --plugin flags, including only plugins whose data files are present.

    These plugins read their own data files (not the VEP cache), so they work in --gtf mode. A
    plugin whose configured path is a placeholder / missing is skipped rather than failing the
    run. SIFT/PolyPhen — and thus Condel — are cache-only and unavailable here; missense is
    covered by AlphaMissense + PrimateAI + CADD.
    """
    pd = VEP_PD
    specs = []
    if _have(pd["cadd_snv"], pd["cadd_indel"]):
        specs.append(f"--plugin CADD,{pd['cadd_snv']},{pd['cadd_indel']}")
    if _have(pd["spliceai_snv"], pd["spliceai_indel"]):
        specs.append(f"--plugin SpliceAI,snv={pd['spliceai_snv']},indel={pd['spliceai_indel']}")
    if _have(pd["primateai"]):
        specs.append(f"--plugin PrimateAI,{pd['primateai']}")
    if _have(pd["alphamissense"]):
        specs.append(f"--plugin AlphaMissense,file={pd['alphamissense']}")
    if _have(pd["loftee_human_ancestor"], pd["loftee_conservation"], pd["loftee_gerp"]):
        specs.append(
            f"--plugin LoF,loftee_path:{pd['loftee_path']},"
            f"human_ancestor_fa:{pd['loftee_human_ancestor']},"
            f"conservation_file:{pd['loftee_conservation']},"
            f"gerp_bigwig:{pd['loftee_gerp']}")
    skipped = {"CADD", "SpliceAI", "PrimateAI", "AlphaMissense", "LoF"} - {
        s.split(",")[0].split()[-1] for s in specs}
    if skipped:
        sys.stderr.write(f"variant_annotation: skipping VEP plugins (no data): {sorted(skipped)}\n")
    return " ".join(specs)


PLUGIN_ARGS = _plugin_args()


rule prepare_gtf:
    """Sort + bgzip + tabix the gencode GTF for VEP `--gtf` (custom annotation)."""
    input:
        gtf=config["gtf"],
    output:
        gtf=REF / "annotation.sorted.gtf.gz",
        tbi=REF / "annotation.sorted.gtf.gz.tbi",
    conda:
        "../../envs/vep.yaml"
    shell:
        r"(zcat {input.gtf} | grep ^'#'; "
        r" zcat {input.gtf} | grep -v ^'#' | sort -k1,1 -k4,4n) "
        r" | bgzip > {output.gtf} && tabix -p gff {output.gtf}"


checkpoint chunk_vcf:
    """Split the fastVEP-filtered VCF into fixed-size bgzipped chunks of
    `additional.vep.variants_per_chunk` variants each (header re-attached, tabix-indexed)."""
    input:
        vcf=FASTVEP_FILTERED_VCF,
    output:
        chunkdir=directory(CHUNKS),
    params:
        region=_REGION,
        n=CHUNK_SIZE,
    conda:
        "../../envs/vep.yaml"
    shell:
        r"""
        mkdir -p {output.chunkdir}
        hdr={output.chunkdir}/.header.vcf
        bcftools view -h {input.vcf} > "$hdr"
        bcftools view -H {params.region} {input.vcf} \
          | split -l {params.n} -d -a 5 - {output.chunkdir}/body_
        for b in {output.chunkdir}/body_*; do
            idx=${{b##*body_}}
            cat "$hdr" "$b" | bgzip > {output.chunkdir}/chunk_$idx.vcf.gz
            tabix -p vcf {output.chunkdir}/chunk_$idx.vcf.gz
            rm -f "$b"
        done
        rm -f "$hdr"
        """


rule vep:
    """VEP + plugins on one variant chunk → tab output. Uploaded_variation =
    chrom_start_ref_alt (join key). Uses --gtf (gencode v34 namespace), not a cache."""
    input:
        vcf=CHUNKS / "chunk_{chunk}.vcf.gz",
        fasta=config["fasta"],
        gtf=REF / "annotation.sorted.gtf.gz",
    output:
        tab=OUT / "vep" / "chunk_{chunk}.vep.tsv.gz",
    params:
        plugin_dir=VEP_CFG["plugin_dir"],
        assembly=config["assembly"],
        fork=VEP_CFG["fork"],
        plugins=PLUGIN_ARGS,
        loftee_path=VEP_PD["loftee_path"],
        distance=VEP_DISTANCE,
    conda:
        "../../envs/vep.yaml"
    shell:
        # PERL5LIB must include LOFTEE for its modules to load; ${{PERL5LIB:-}} guards `set -u`.
        # NB: NO --offline (forces a cache, incompatible with --gtf). SIFT/PolyPhen + cached
        # gnomAD AF are cache-only → unavailable in --gtf mode.
        "PERL5LIB={params.loftee_path}:${{PERL5LIB:-}} "
        "vep --gtf {input.gtf} --fasta {input.fasta} "
        "    --dir_plugins {params.plugin_dir} --assembly {params.assembly} "
        "    --fork {params.fork} --force_overwrite --tab --compress_output bgzip --no_stats "
        "    --symbol --canonical --biotype --numbers --distance {params.distance} "
        "    --input_file {input.vcf} --output_file {output.tab} "
        "    {params.plugins}"


rule parse_vep:
    """VEP tab → tidy per-chunk parquet (variant×transcript), carrying variant_id + transcript."""
    input:
        tab=OUT / "vep" / "chunk_{chunk}.vep.tsv.gz",
    output:
        parquet=OUT / "vep" / "chunk_{chunk}.vep.parquet",
    conda:
        "../../envs/parse.yaml"
    script:
        "../scripts/parse_vep.py"


def _vep_chunk_parquets(wildcards):
    return expand(str(OUT / "vep" / "chunk_{chunk}.vep.parquet"), chunk=chunk_ids())


rule vep_parquet:
    """Concat the per-chunk VEP parquets -> parts/vep.parquet (flat; combine joins on variant_id+transcript)."""
    input:
        parquets=_vep_chunk_parquets,
    output:
        parquet=PARTS / "vep.parquet",
    run:
        import polars as pl
        (pl.concat([pl.scan_parquet(p) for p in input.parquets])
           .collect(engine="streaming")
           .write_parquet(output.parquet))
