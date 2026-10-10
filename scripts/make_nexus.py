#!/usr/bin/env python3
"""
Write one NEXUS file holding a protein alignment (TAXA + CHARACTERS blocks) and a tree (TREES block),
e.g. for viewing alignment and tree side by side. Taxa are the intersection of alignment sequences
and tree tips; the tree is pruned to them if needed. Support values in the tree are kept.

Usage:
    make_nexus.py <alignment.fa> <tree.newick> <out.nex>
"""
import sys
from io import StringIO

from Bio import Phylo, SeqIO


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    aln_path, tree_path, out_path = sys.argv[1:]

    aln = {r.id: str(r.seq).upper() for r in SeqIO.parse(aln_path, "fasta")}
    tree = Phylo.read(tree_path, "newick")
    tips = [t.name for t in tree.get_terminals()]
    common = [t for t in tips if t in aln]
    if not common:
        sys.exit("no taxa shared between alignment and tree")
    for t in tips:
        if t not in aln:
            tree.prune(t)
    lengths = {len(aln[t]) for t in common}
    if len(lengths) != 1:
        sys.exit("alignment sequences differ in length")
    nchar = lengths.pop()

    buf = StringIO()
    Phylo.write(tree, buf, "newick")
    newick = buf.getvalue().strip()

    w = max(len(t) for t in common) + 2
    with open(out_path, "w") as out:
        out.write("#NEXUS\n\n")
        out.write(f"BEGIN TAXA;\n\tDIMENSIONS NTAX={len(common)};\n\tTAXLABELS\n")
        out.writelines(f"\t\t{t}\n" for t in common)
        out.write("\t;\nEND;\n\n")
        out.write(f"BEGIN CHARACTERS;\n\tDIMENSIONS NCHAR={nchar};\n"
                  "\tFORMAT DATATYPE=PROTEIN MISSING=? GAP=-;\n\tMATRIX\n")
        out.writelines(f"\t{t.ljust(w)}{aln[t]}\n" for t in common)
        out.write("\t;\nEND;\n\n")
        out.write(f"BEGIN TREES;\n\tTREE tree1 = {newick}\nEND;\n")
    print(f"[make_nexus] {len(common)} taxa, {nchar} columns -> {out_path}"
          + (f" ({len(tips) - len(common)} tree tips without sequence pruned)" if len(tips) > len(common) else ""),
          file=sys.stderr)


if __name__ == "__main__":
    main()
