#!/bin/bash
# runjob.sh <tag> <seed> <beam nu|nubar> <program args...>  (event mode: slims output; total modes: copies txt)
# Paths are relative to validation/asymmetric-nuclei/. Optional env: BIN (generator binary), Q2P2H_OVERRIDE (MeV),
# DATA_DIR (npz/txt output, default data/), RUNS_DIR (per-job work dirs and logs, default runs/prod).
V="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
tag=$1; seed=$2; beam=$3; shift 3
DATA_DIR=${DATA_DIR:-$V/data}; RUNS_DIR=${RUNS_DIR:-$V/runs/prod}
d=$RUNS_DIR/$tag; mkdir -p "$d" "$DATA_DIR"; cd "$d" || exit 1
echo "seed=$seed args=$* Q2P2H_OVERRIDE=${Q2P2H_OVERRIDE:-none}" > cmd.txt
start=$(date +%s)
MC_SEED=$seed "${BIN:-$V/generator/build/program}" "$@" > out.log 2> err.log; rc=$?
echo "rc=$rc secs=$(( $(date +%s)-start ))" >> cmd.txt
shopt -s nullglob; ev=(events*.root); tot=(total*.txt)
if (( ${#ev[@]} )); then python3 "$V/scripts/generation/slim.py" "${ev[@]}" "$DATA_DIR/$tag.npz" "$beam" >> cmd.txt 2>&1 && rm -f "${ev[@]}"; fi
if (( ${#tot[@]} )); then cp "${tot[@]}" "$DATA_DIR/$tag.txt"; fi
echo "$tag rc=$rc" >> "$RUNS_DIR/done.txt"
