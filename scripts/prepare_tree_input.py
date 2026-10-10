#!/usr/bin/env python3
"""
Pool candidate protein sets (e.g. sorted/<gene>/synth_ready/*.prot.fa from one or more runs) and the
reference proteins into tree-ready FASTA files with names that are safe for IQ-TREE, newick and NEXUS
([A-Za-z0-9_.] only).

Candidate headers written by this pipeline (species|accession|chrom:start-end(strand)|clade|tier)
become <TAG>_<species>[_<accession>]_<chrom>_<start>_<end>_<plus|minus>; any other header is cleaned
as a whole and prefixed with <TAG>_. Reference names get the prefix REF_. Terminal '*' is removed.

A candidate whose protein AND genomic start-end are identical to one already kept (the same
assembly listed in two sets) is skipped and logged in name_map.tsv.

--outgroup FILE adds outgroup proteins (e.g. a paralogue) named OUT_<id>; they are aligned like the
candidates and their tree names are listed in outgroup_names.txt (used by build_tree.sh to root).

Usage:
    prepare_tree_input.py --reference_proteins REF.faa --set TAG=FILE [--set TAG=FILE ...] \
                          [--outgroup FILE] -o OUTDIR

Writes OUTDIR/refs.prot.fa, OUTDIR/candidates.prot.fa, OUTDIR/name_map.tsv,
OUTDIR/outgroup_names.txt
"""
import argparse
import os
import re
import sys

from Bio import SeqIO


def clean(s):
    s = s.replace("(+)", "_plus").replace("(-)", "_minus")
    return re.sub(r"_+", "_", re.sub(r"[^A-Za-z0-9.]", "_", s)).strip("_")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference_proteins", required=True)
    ap.add_argument("--set", action="append", required=True, metavar="TAG=FILE",
                    help="candidate protein FASTA with a short tag used as name prefix (repeatable)")
    ap.add_argument("--outgroup", help="outgroup protein FASTA (named OUT_<id>, used to root)")
    ap.add_argument("-o", "--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    seen, kept = set(), {}
    n_ref = n_cand = n_dup = 0
    with open(os.path.join(args.outdir, "name_map.tsv"), "w") as m, \
         open(os.path.join(args.outdir, "refs.prot.fa"), "w") as rf, \
         open(os.path.join(args.outdir, "candidates.prot.fa"), "w") as cf:
        m.write("tree_name\tset\toriginal_header\n")
        for rec in SeqIO.parse(args.reference_proteins, "fasta"):
            name = "REF_" + clean(rec.id)
            if name in seen:
                sys.exit(f"duplicate reference name after cleaning: {name}")
            seen.add(name)
            rf.write(f">{name}\n{str(rec.seq).rstrip('*')}\n")
            m.write(f"{name}\treference\t{rec.id}\n")
            n_ref += 1
        for spec in args.set:
            tag, path = spec.split("=", 1)
            for rec in SeqIO.parse(path, "fasta"):
                seq = str(rec.seq).rstrip("*")
                parts = rec.id.split("|")
                coords = None
                if len(parts) >= 3 and ":" in parts[2]:
                    species, acc, loc = parts[:3]
                    coords = loc.split(":", 1)[1]
                    core = species if species == acc else f"{species}_{acc}"
                    name = f"{tag}_{clean(core)}_{clean(loc)}"
                else:
                    name = f"{tag}_{clean(rec.id)}"
                if coords is not None and (seq, coords) in kept:
                    m.write(f"(skipped: duplicate of {kept[(seq, coords)]})\t{tag}\t{rec.id}\n")
                    n_dup += 1
                    continue
                base, k = name, 2
                while name in seen:
                    name, k = f"{base}_{k}", k + 1
                seen.add(name)
                if coords is not None:
                    kept[(seq, coords)] = name
                cf.write(f">{name}\n{seq}\n")
                m.write(f"{name}\t{tag}\t{rec.id}\n")
                n_cand += 1
        out_names = []
        if args.outgroup:
            for rec in SeqIO.parse(args.outgroup, "fasta"):
                name = "OUT_" + clean(rec.id)
                if name in seen:
                    sys.exit(f"duplicate outgroup name after cleaning: {name}")
                seen.add(name)
                cf.write(f">{name}\n{str(rec.seq).rstrip('*')}\n")
                m.write(f"{name}\toutgroup\t{rec.id}\n")
                out_names.append(name)
    with open(os.path.join(args.outdir, "outgroup_names.txt"), "w") as fh:
        fh.writelines(f"{n}\n" for n in out_names)
    print(f"[prepare_tree_input] {n_ref} references, {n_cand} candidates, {len(out_names)} outgroups, "
          f"{n_dup} duplicates skipped",
          file=sys.stderr)


if __name__ == "__main__":
    main()
