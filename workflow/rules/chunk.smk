"""Shared VCF chunking — split fastVEP's filtered VCF into fixed-size chunks.

Used by BOTH VEP and NMD (deeprvat-style scatter: even parallelism, bounded per-job memory — a
single job over the whole 20M-variant set OOMs/straggles). Included when either additional.vep or
additional.nmd is enabled. Chunk size = additional.variants_per_chunk (shared; default 500k).
"""
import os

CHUNKS = OUT / "chunks"                  # chunk_vcf writes chunk_NNNNN.vcf.gz here (temp; VEP + NMD consume)
CHUNK_SIZE = int(config["additional"].get("variants_per_chunk", 500_000))
# Optional chromosome filter applied before chunking ([] / unset = the whole filtered VCF).
_REGION = ("-r " + ",".join(config["chromosomes"])) if config.get("chromosomes") else ""


def chunk_ids():
    """Resolve the dynamic chunk wildcards produced by the chunk_vcf checkpoint."""
    cdir = checkpoints.chunk_vcf.get().output.chunkdir
    return sorted(glob_wildcards(os.path.join(cdir, "chunk_{chunk}.vcf.gz")).chunk)


checkpoint chunk_vcf:
    """Split fastvep/small.vcf.gz into fixed-size bgzipped chunks (header re-attached, tabix'd).

    Uses envs/fastvep.yaml (bcftools + htslib) — always built since fastVEP runs — so enabling NMD
    without VEP doesn't force building the heavy ensembl-vep env just to split the VCF."""
    input:
        vcf=FASTVEP_FILTERED_VCF,
    output:
        chunkdir=temp(directory(CHUNKS)),   # intermediate — deleted once VEP + NMD have consumed the chunks
    params:
        region=_REGION,
        n=CHUNK_SIZE,
    conda:
        "../../envs/fastvep.yaml"
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
