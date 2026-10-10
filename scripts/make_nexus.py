#!/usr/bin/env python3
"""
Write one NEXUS file holding a protein alignment (TAXA + CHARACTERS blocks) and a tree (TREES block),
e.g. for viewing alignment and tree side by side. Taxa are the intersection of alignment sequences
and tree tips; the tree is pruned to them if needed. Support values in the tree are kept.

With --outgroup NAMES_FILE (one tip name per line, e.g. outgroup_names.txt from
prepare_tree_input.py) the tree is rooted on the branch separating the outgroup from all other tips,
giving a bifurcating root (IQ-TREE -o only draws the outgroup at a basal multifurcation). Support
values stay attached to the same bipartitions; this is checked and the script fails otherwise.

Usage:
    make_nexus.py <alignment.fa> <tree.newick> <out.nex> [--outgroup NAMES_FILE]
"""
import argparse
import sys
from io import StringIO

from Bio import Phylo, SeqIO


def bipartition_support(tree):
    """{frozenset(smaller side of the split): support} for every internal branch with a support value"""
    allt = frozenset(t.name for t in tree.get_terminals())
    out = {}
    for c in tree.get_nonterminals():
        if c is tree.root or c.confidence is None:
            continue
        side = frozenset(t.name for t in c.get_terminals())
        other = allt - side
        key = min(side, other, key=lambda x: (len(x), sorted(x)))
        out[key] = c.confidence
    return out


def root_on_outgroup(tree, outgroup):
    tips = [t.name for t in tree.get_terminals()]
    out = [t for t in tips if t in outgroup]
    ingroup = [t for t in tips if t not in outgroup]
    if not out or not ingroup:
        sys.exit("outgroup: none of the outgroup names are in the tree")
    before = bipartition_support(tree)
    # move the root into the outgroup first, then onto the stem of the ingroup
    tree.root_with_outgroup(out[0])
    ca = tree.common_ancestor(ingroup)
    if set(t.name for t in ca.get_terminals()) != set(ingroup):
        sys.exit("outgroup: ingroup is not monophyletic with respect to the outgroup; not rooting")
    tree.root_with_outgroup(ca)   # the ingroup MRCA is now the (multifurcating) root node
    root = tree.root
    out_kids = [c for c in root.clades if {t.name for t in c.get_terminals()} <= set(out)]
    in_kids = [c for c in root.clades if c not in out_kids]
    if not out_kids or not in_kids:
        sys.exit("outgroup: could not separate outgroup and ingroup at the root")
    if len(out_kids) == 1:
        out_node = out_kids[0]
    else:
        out_node = Phylo.BaseTree.Clade(branch_length=0.0, clades=out_kids)
    stem = out_node.branch_length or 0.0
    out_node.branch_length = stem / 2
    in_node = Phylo.BaseTree.Clade(branch_length=stem / 2, clades=in_kids, confidence=out_node.confidence)
    tree.root = Phylo.BaseTree.Clade(clades=[in_node, out_node])
    tree.rooted = True
    after = bipartition_support(tree)
    lost = {k: v for k, v in before.items() if k in after and after[k] != v}
    if lost:
        sys.exit(f"outgroup: support values moved during rerooting for {len(lost)} branches")
    print(f"[make_nexus] rooted between {len(out)} outgroup and {len(ingroup)} ingroup tips; "
          f"{sum(k in after for k in before)}/{len(before)} supported branches checked", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("alignment")
    ap.add_argument("tree")
    ap.add_argument("out")
    ap.add_argument("--outgroup", help="file with outgroup tip names, one per line")
    args = ap.parse_args()
    aln_path, tree_path, out_path = args.alignment, args.tree, args.out

    aln = {r.id: str(r.seq).upper() for r in SeqIO.parse(aln_path, "fasta")}
    tree = Phylo.read(tree_path, "newick")
    tips = [t.name for t in tree.get_terminals()]
    common = [t for t in tips if t in aln]
    if not common:
        sys.exit("no taxa shared between alignment and tree")
    for t in tips:
        if t not in aln:
            tree.prune(t)
    if args.outgroup:
        names = {l.strip() for l in open(args.outgroup) if l.strip()}
        if names:
            root_on_outgroup(tree, names)
            common = [t.name for t in tree.get_terminals()]
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
        rooted = "[&R] " if args.outgroup else ""
        out.write(f"BEGIN TREES;\n\tTREE tree1 = {rooted}{newick}\nEND;\n")
    print(f"[make_nexus] {len(common)} taxa, {nchar} columns -> {out_path}"
          + (f" ({len(tips) - len(common)} tree tips without sequence pruned)" if len(tips) > len(common) else ""),
          file=sys.stderr)


if __name__ == "__main__":
    main()
