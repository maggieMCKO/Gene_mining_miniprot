# Changelog (fork)

Changes relative to MaximePolicarpo/Gene_mining_miniprot @ 2ec3b67.

## Unreleased — tuned-tiers

### Added
- `mine_gene.sh --miniprot-args STR`: extra miniprot arguments for both passes (recommended `-B 15`).
- `mine_gene.sh --dedup-prefer-intact` / `dedup_loci.py --prefer-intact`: per locus, prefer an intact
  model (no internal stop / frameshift) that is ≥ 90 % as long as the longest model there.
- `scripts/assess_tiers.py`: per-locus `tier` (`synth_ready`, `synth_ready_extended`, `asr_usable`,
  `not_usable`) and its criteria columns appended to `summary.tsv`; writes
  `<prefix>.tiers.{cds,prot}.fa`. In-frame start/stop extension (≤ 10 codons, never through N).
- `scripts/sort_sequences.py`: pools per-genome FASTAs into `sorted/<gene|suspect>/<call or tier>/`;
  called at the end of `run_batch.sh`.

### Unchanged
- Default behaviour of every original option and step.
