"""Drive an AbSplice-DNA run and aggregate to a per-variant parquet.

This is the integration point with gagneurlab/absplice. AbSplice has its own snakemake
workflow that needs SpliceMap reference data per tissue; rather than vendoring it, point
`--absplice-config` at your AbSplice checkout's config and let it run, then aggregate its
tissue outputs here.

Expected output parquet schema (keyed on the variant):
    chrom, start, ref, alt, AbSplice_DNA            # max across requested tissues
    [optional] AbSplice_DNA_<tissue> per tissue     # wide, if you want per-tissue

Implement the body against your AbSplice install. Kept as an explicit stub so the DAG is
wired but the heavy AbSplice setup is opt-in.
"""
from __future__ import annotations

import argparse
import sys


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--absplice-config", required=True)
    ap.add_argument("--tissues", required=True, help="comma-separated tissue names")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    sys.exit(
        "run_absplice.py is a stub — wire it to your gagneurlab/absplice checkout.\n"
        f"  vcf={args.vcf}\n  absplice_config={args.absplice_config}\n"
        f"  tissues={args.tissues}\n  out={args.out}\n"
        "Run AbSplice on the VCF, then aggregate its per-tissue AbSplice_DNA scores to a\n"
        "parquet keyed on (chrom,start,ref,alt) with an AbSplice_DNA column."
    )


if __name__ == "__main__":
    main()
