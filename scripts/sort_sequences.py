#!/usr/bin/env python3
"""
Collect per-genome mine_gene.sh output sequences into folders by functional call
and tree placement:

    <outdir>/<GENE>/<functional_call>/    clade_call == GENE
    <outdir>/suspect/<functional_call>/   OUTGROUP_SUSPECT / AMBIGUOUS(...) / NOT_IN_TREE

Each folder gets <GENE>.<group>.<call>.cds.fa, .prot.fa and .summary.tsv
(the matching rows of the per-genome summary.tsv files).

If the per-genome outputs have usability tiers (*.tiers.{cds,prot}.fa, written by
scripts/assess_tiers.py), three more folders are written per group:

    <group>/synth_ready/            tier synth_ready (CDS as predicted, ATG..stop)
    <group>/synth_ready_extended/   tier synth_ready_extended (CDS incl. the in-frame extension)
    <group>/asr_usable/             all in-frame models: synth_ready + synth_ready_extended +
                                    asr_usable; terminal stop codon removed (for HyPhy/FastML)

Usage:
    sort_sequences.py <results_dir> <gene_name> [--calls functional_complete pseudogene] [--outdir DIR]
"""
import argparse
import glob
import os
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results_dir")
    ap.add_argument("gene_name")
    ap.add_argument("--calls", nargs="+", default=["functional_complete", "pseudogene"])
    ap.add_argument("--outdir", default=None, help="default: <results_dir>/sorted")
    args = ap.parse_args()

    gene = args.gene_name
    outdir = args.outdir or os.path.join(args.results_dir, "sorted")
    gene_dir = os.path.join(args.results_dir, gene)

    def group_of(clade_call):
        return gene if clade_call == gene else "suspect"

    # summary rows, keyed by gene_id (= FASTA header minus the |clade|call suffix)
    header_cols, rows = None, {}
    for f in sorted(glob.glob(os.path.join(gene_dir, "*", f"*.{gene}.summary.tsv"))):
        with open(f) as fh:
            first = fh.readline().rstrip("\n")
            if not first:
                continue
            header_cols = first.split("\t")
            for line in fh:
                vals = line.rstrip("\n").split("\t")
                rows[vals[0]] = vals

    buckets = {}  # (call, group) -> {"cds": [...], "prot": [...], "rows": [...]}
    for kind in ("cds", "prot"):
        pattern = os.path.join(gene_dir, "*", f"*.{gene}*.{kind}.fa")
        for f in sorted(glob.glob(pattern)):
            if f.endswith(f".tiers.{kind}.fa"):
                continue
            for rec in SeqIO.parse(f, "fasta"):
                parts = rec.id.split("|")
                gene_id, clade_call, call = "|".join(parts[:-2]), parts[-2], parts[-1]
                if call not in args.calls:
                    continue
                b = buckets.setdefault((call, group_of(clade_call)), {"cds": [], "prot": [], "rows": []})
                b[kind].append(rec)
                if kind == "cds":
                    b["rows"].append(rows[gene_id])

    # usability tiers
    tier_dirs = {"synth_ready": ["synth_ready"],
                 "synth_ready_extended": ["synth_ready_extended"],
                 "asr_usable": ["asr_usable", "synth_ready", "synth_ready_extended"]}
    have_tiers = False
    for kind in ("cds", "prot"):
        for f in sorted(glob.glob(os.path.join(gene_dir, "*", f"*.{gene}.tiers.{kind}.fa"))):
            have_tiers = True
            for rec in SeqIO.parse(f, "fasta"):
                parts = rec.id.split("|")
                gene_id, clade_call, tier = "|".join(parts[:-2]), parts[-2], parts[-1]
                for folder, members in tier_dirs.items():
                    if tier not in members:
                        continue
                    r = rec
                    if folder == "asr_usable":
                        seq = str(rec.seq)
                        if kind == "cds" and seq[-3:].upper() in ("TAA", "TAG", "TGA"):
                            seq = seq[:-3]
                        if kind == "prot":
                            seq = seq.rstrip("*")
                        r = SeqRecord(Seq(seq), id=rec.id, description="")
                    b = buckets.setdefault((folder, group_of(clade_call)), {"cds": [], "prot": [], "rows": []})
                    b[kind].append(r)
                    if kind == "cds":
                        b["rows"].append(rows[gene_id])
    calls_out = list(args.calls) + (list(tier_dirs) if have_tiers else [])

    for group in (gene, "suspect"):
        for call in calls_out:
            b = buckets.get((call, group), {"cds": [], "prot": [], "rows": []})
            d = os.path.join(outdir, group, call)
            os.makedirs(d, exist_ok=True)
            prefix = os.path.join(d, f"{gene}.{group}.{call}")
            for kind in ("cds", "prot"):
                SeqIO.write(b[kind], f"{prefix}.{kind}.fa", "fasta-2line")
            with open(f"{prefix}.summary.tsv", "w") as out:
                out.write("\t".join(header_cols) + "\n")
                for r in b["rows"]:
                    out.write("\t".join(r) + "\n")
            print(f"[sort] {group}/{call}: {len(b['cds'])} cds, {len(b['prot'])} prot", file=sys.stderr)


if __name__ == "__main__":
    main()
