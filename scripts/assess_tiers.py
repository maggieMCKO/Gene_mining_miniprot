#!/usr/bin/env python3
"""
Assign each final locus a usability tier, append the tier columns to <prefix>.summary.tsv,
and write the final per-locus sequences to <prefix>.tiers.{cds,prot}.fa.

Tiers (one per locus):
  synth_ready           full-length CDS fit for synthesis, as predicted by miniprot: ATG start,
                        terminal stop, no internal stop / frameshift / in-frame stop, no N or other
                        ambiguous base anywhere in the CDS, both
                        reference termini covered (<= --end-tol aa missing in the reference core),
                        no internal gap >= --gap-len aa that no reference shares, all introns
                        canonical (GT/GC...AG) and >= --min-intron bp
  synth_ready_extended  as synth_ready, but the start ATG and/or the stop codon was found only by
                        extending the CDS in frame along the genome, by at most --ext-codons codons,
                        never through N/ambiguous bases (codon counts, incl. the stop codon, in
                        ext_codons_N / ext_codons_C). The written CDS includes the
                        extension.
  asr_usable            in-frame model usable in a codon alignment (aBSREL, FastML, ...): no internal
                        stop / frameshift / in-frame stop, >= --min-cov of the reference core covered.
                        Partial models, gaps and Ns are allowed.
  not_usable            everything else

Usage:
    assess_tiers.py <spdir> <prefix> <genome.fa> <reference_proteins.fa> [options]
"""
import argparse
import csv
import os
import re

from Bio import SeqIO
from Bio.Seq import Seq

STOPS = {"TAA", "TAG", "TGA"}
NEW_COLS = ["start_M", "end_stop", "ext_codons_N", "ext_codons_C", "n_X", "n_ambiguous_nt", "ref_core_missing_N",
            "ref_core_missing_C", "max_unexplained_gap", "ref_core_cov", "min_intron",
            "noncanonical_splices", "tier"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("spdir")
    ap.add_argument("prefix", help="output prefix, e.g. <spdir>/<label>.<gene>")
    ap.add_argument("genome")
    ap.add_argument("reference_proteins")
    ap.add_argument("--ext-codons", type=int, default=10)
    ap.add_argument("--end-tol", type=int, default=5)
    ap.add_argument("--gap-len", type=int, default=20)
    ap.add_argument("--min-cov", type=float, default=0.5)
    ap.add_argument("--min-intron", type=int, default=60)
    args = ap.parse_args()

    tmp = os.path.join(args.spdir, "tmp")
    summ_path = f"{args.prefix}.summary.tsv"
    with open(summ_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames
        summary = list(reader)
    cds_out = open(f"{args.prefix}.tiers.cds.fa", "w")
    prot_out = open(f"{args.prefix}.tiers.prot.fa", "w")
    if not summary:
        cds_out.close(); prot_out.close()
        return

    ref_ids = [r.id for r in SeqIO.parse(args.reference_proteins, "fasta")]
    cds = {r.id: str(r.seq).upper() for r in SeqIO.parse(os.path.join(tmp, "dedup.cds.fa"), "fasta")}
    aln = {r.id: str(r.seq) for r in SeqIO.parse(os.path.join(tmp, "placed.aln"), "fasta")}
    want = {s["mrna_id"] for s in summary}
    blocks = {}
    for line in open(os.path.join(tmp, "raw.gff")):
        c = line.rstrip("\n").split("\t")
        if len(c) > 8 and c[2] == "CDS":
            m = re.search(r"Parent=([^;]+)", c[8])
            if m and m.group(1) in want:
                blocks.setdefault(m.group(1), []).append((c[0], int(c[3]), int(c[4]), c[6]))
    genome = SeqIO.index(args.genome, "fasta")

    # reference core = alignment columns where >= 4 references (or all, if fewer) have a residue
    refs_in = [r for r in ref_ids if r in aln]
    need = min(4, len(refs_in))
    L = len(aln[refs_in[0]])
    core = [i for i in range(L) if sum(aln[r][i] != "-" for r in refs_in) >= need]
    ref_core = ["".join(aln[r][i] for i in core) for r in refs_in]

    def gap_explained(a, b):
        return any(rc[a:b].count("-") >= 0.8 * (b - a) for rc in ref_core)

    for s in summary:
        mid = s["mrna_id"]
        seq = cds[mid]
        bl = sorted(blocks.get(mid, []), key=lambda b: b[1])
        chrom, strand = bl[0][0], bl[0][3]
        g = genome[chrom].seq

        prot = str(Seq(seq[: len(seq) // 3 * 3]).translate())
        start_M, end_stop = prot[:1] == "M", prot.endswith("*")

        # in-frame extension (only used for the synth_ready_extended tier)
        n = args.ext_codons * 3
        first, last = bl[0][1], bl[-1][2]
        if strand == "+":
            down = str(g[last: last + n]).upper()
            up = str(g[max(0, first - 1 - n): first - 1]).upper()
        else:
            down = str(g[max(0, first - 1 - n): first - 1].reverse_complement()).upper()
            up = str(g[last: last + n].reverse_complement()).upper()
        ext_c, ext_n, ok_c, ok_n = "", "", end_stop, start_M
        if not end_stop:
            for i in range(0, len(down) - 2, 3):
                ext_c += down[i:i + 3]
                if down[i:i + 3] in STOPS:
                    ok_c = True
                    break
            if not ok_c:
                ext_c = ""
        if not start_M:
            up = up[len(up) % 3:]
            for i in range(len(up) - 3, -1, -3):
                codon = up[i:i + 3]
                if codon in STOPS:
                    break
                if codon == "ATG":
                    ext_n, ok_n = up[i:], True
                    break
        if set(ext_c) - set("ACGT"):   # never extend through assembly gaps / ambiguous bases
            ext_c, ok_c = "", end_stop
        if set(ext_n) - set("ACGT"):
            ext_n, ok_n = "", start_M
        ext_seq = ext_n + seq + ext_c

        # coverage and gaps against the reference core
        cand = "".join(aln[mid][i] for i in core) if mid in aln else "-" * len(core)
        lead = len(cand) - len(cand.lstrip("-"))
        trail = len(cand) - len(cand.rstrip("-"))
        core_cov = 1 - cand.count("-") / len(core)
        max_gap = 0
        for m in re.finditer(r"-+", cand):
            a, b = m.start(), m.end()
            if a > 0 and b < len(cand) and not gap_explained(a, b):
                max_gap = max(max_gap, b - a)

        # intron sanity
        introns = [(bl[j][2] + 1, bl[j + 1][1] - 1) for j in range(len(bl) - 1)]
        min_intron = min((e - b0 + 1 for b0, e in introns), default=0)
        noncanon = 0
        for b0, e in introns:
            intr = g[b0 - 1:e]
            if strand == "-":
                intr = intr.reverse_complement()
            intr = str(intr).upper()
            if not (intr[:2] in ("GT", "GC") and intr[-2:] == "AG"):
                noncanon += 1

        n_X = prot.count("X")
        n_amb = sum(ch not in "ACGT" for ch in seq)   # an N can hide in a synonymous codon (CTN = L)
        intact = int(s["internal_stops"]) == 0 and int(s["gff_frameshift"]) == 0 \
            and int(s["gff_stopcodon"]) == 0
        body_ok = intact and n_X == 0 and n_amb == 0 and lead <= args.end_tol and trail <= args.end_tol \
            and max_gap < args.gap_len and len(seq) % 3 == 0 and noncanon == 0 \
            and (not introns or min_intron >= args.min_intron)
        if body_ok and start_M and end_stop:
            tier, out_seq = "synth_ready", seq
        elif body_ok and ok_n and ok_c:
            tier, out_seq = "synth_ready_extended", ext_seq
        elif intact and core_cov >= args.min_cov:
            tier, out_seq = "asr_usable", seq
        else:
            tier, out_seq = "not_usable", seq

        s.update({"start_M": start_M, "end_stop": end_stop,
                  "ext_codons_N": len(ext_n) // 3 if tier == "synth_ready_extended" else 0,
                  "ext_codons_C": len(ext_c) // 3 if tier == "synth_ready_extended" else 0,
                  "n_X": n_X, "n_ambiguous_nt": n_amb, "ref_core_missing_N": lead, "ref_core_missing_C": trail,
                  "max_unexplained_gap": max_gap, "ref_core_cov": round(core_cov, 3),
                  "min_intron": min_intron, "noncanonical_splices": noncanon, "tier": tier})
        hdr = f">{s['gene_id']}|{s['clade_call']}|{tier}"
        cds_out.write(f"{hdr}\n{out_seq}\n")
        prot_out.write(f"{hdr}\n{str(Seq(out_seq[: len(out_seq) // 3 * 3]).translate())}\n")

    cds_out.close(); prot_out.close()
    cols = [c for c in header if c not in NEW_COLS] + NEW_COLS
    with open(summ_path, "w") as out:
        out.write("\t".join(cols) + "\n")
        for s in summary:
            out.write("\t".join(str(s[c]) for c in cols) + "\n")


if __name__ == "__main__":
    main()
