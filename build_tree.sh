#!/usr/bin/env bash
# build_tree.sh - protein alignment + ML tree + combined NEXUS for candidate sets produced by
# mine_gene.sh / run_batch.sh (e.g. sorted/<gene>/synth_ready/*.prot.fa), together with the
# reference proteins.
#
#   1. prepare_tree_input.py: pool sets, tree-safe names (name_map.tsv), drop cross-set duplicates
#   2. muscle -align the reference proteins (the alignment frame)
#   3. mafft --add candidates --keeplength onto that frame
#   4. trimal -gt 0.7 -cons 60
#   5. iqtree -m TEST -bb 1000 (ModelFinder + 1000 ultrafast bootstraps)
#   6. make_nexus.py: trimmed alignment + tree in one NEXUS file
#
# Usage:
#   build_tree.sh --reference_proteins REF.faa --set TAG=PROT.fa [--set TAG=PROT.fa ...] \
#                 -o OUTDIR [--prefix NAME] [-t THREADS] [--trimal-args "..."] [--iqtree-args "..."]
#
# Required:
#   --reference_proteins FILE   unaligned reference proteins (e.g. the mine_gene.sh ingroup set)
#   --set TAG=FILE              candidate proteins; TAG (short, e.g. NEC) prefixes the tree names
#   -o/--outdir DIR
# Optional:
#   --prefix NAME               output file prefix [tree]
#   -t/--threads INT            [4]
#   --trimal-args STR           [-gt 0.7 -cons 60]
#   --iqtree-args STR           [-m TEST -bb 1000]
#
# Note: --keeplength drops residues that would be insertions relative to the reference frame, and
# the tree is unrooted unless you add outgroup sequences as one more --set.
#
# Environment: environment_tree.yml (muscle, mafft, trimal, iqtree, biopython).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF=""; OUTDIR=""; PREFIX="tree"; THREADS=4
TRIMAL_ARGS="-gt 0.7 -cons 60"; IQTREE_ARGS="-m TEST -bb 1000"; SETS=()

usage() { grep '^# ' "$0" | sed 's/^# \?//'; exit 1; }
while [[ $# -gt 0 ]]; do
    case "$1" in
        --reference_proteins) REF="$2"; shift 2 ;;
        --set) SETS+=(--set "$2"); shift 2 ;;
        -o|--outdir) OUTDIR="$2"; shift 2 ;;
        --prefix) PREFIX="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        --trimal-args) TRIMAL_ARGS="$2"; shift 2 ;;
        --iqtree-args) IQTREE_ARGS="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1" >&2; usage ;;
    esac
done
[[ -z "$REF" || -z "$OUTDIR" || ${#SETS[@]} -eq 0 ]] && usage
for tool in muscle mafft trimal iqtree python3; do
    command -v "$tool" >/dev/null || { echo "ERROR: $tool not on PATH (see environment_tree.yml)" >&2; exit 1; }
done

log() { echo "[$(date +%H:%M:%S)] $*" >&2; }
mkdir -p "$OUTDIR"

log "Preparing input"
python3 "${SCRIPT_DIR}/scripts/prepare_tree_input.py" --reference_proteins "$REF" "${SETS[@]}" -o "$OUTDIR"
cd "$OUTDIR"

log "Aligning reference proteins (muscle)"
muscle -align refs.prot.fa -output refs.prot.aln 2> muscle.log

log "Adding candidates onto the reference frame (mafft --add --keeplength)"
mafft --thread "$THREADS" --add candidates.prot.fa --keeplength refs.prot.aln \
    > "${PREFIX}.prot.aln" 2> mafft.log

log "Trimming (trimal ${TRIMAL_ARGS})"
trimal -in "${PREFIX}.prot.aln" -out "${PREFIX}.prot.aln.trimmed" $TRIMAL_ARGS > trimal.log 2>&1

log "Building tree (iqtree ${IQTREE_ARGS})"
iqtree -s "${PREFIX}.prot.aln.trimmed" --seqtype AA $IQTREE_ARGS -nt "$THREADS" \
    --prefix "$PREFIX" -redo --quiet

log "Writing NEXUS (alignment + tree)"
python3 "${SCRIPT_DIR}/scripts/make_nexus.py" "${PREFIX}.prot.aln.trimmed" "${PREFIX}.treefile" \
    "${PREFIX}.prot.aln_with_tree.nex"

log "Done: ${OUTDIR}/${PREFIX}.prot.aln_with_tree.nex"
