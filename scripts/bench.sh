#!/usr/bin/env bash
# ============================================================
# End-to-end wall-clock benchmarks using hyperfine.
#
# Prerequisites:
#   cargo build --release
#   cargo run --release --bin gen_fastq -- --help
#   hyperfine (https://github.com/sharkdp/hyperfine)
#
# Usage:
#   bash scripts/bench.sh [--reads N] [--read-len L] [--threads T]
#
# Example:
#   bash scripts/bench.sh --reads 10000000 --read-len 150 --threads 8
# ============================================================

set -euo pipefail

READS=5000000
READ_LEN=150
NUM_SEQS=10000
THREADS=""   # empty = use all CPUs

# Parse optional args
while [[ $# -gt 0 ]]; do
  case "$1" in
    --reads)    READS="$2";    shift 2 ;;
    --read-len) READ_LEN="$2"; shift 2 ;;
    --threads)  THREADS="$2";  shift 2 ;;
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
done

BIN="./target/release/mmfqcount"
GEN="./target/release/gen_fastq"
TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

echo "=== Building release binary ==="
cargo build --release --quiet

echo ""
echo "=== Generating test data ($READS reads, len=$READ_LEN, pool=$NUM_SEQS) ==="

R1="$TMPDIR/r1.fastq"
R2="$TMPDIR/r2.fastq"
R1_GZ="$TMPDIR/r1.fastq.gz"
OUT="$TMPDIR/out.tsv"

"$GEN" --reads "$READS" --read-len "$READ_LEN" --num-seqs "$NUM_SEQS" \
       --output "$R1" --output-r2 "$R2"

# Also produce a gzip version for the gz benchmark
"$GEN" --reads "$READS" --read-len "$READ_LEN" --num-seqs "$NUM_SEQS" \
       --output "$R1_GZ"

FILE_SIZE_MB=$(du -sm "$R1" | cut -f1)
echo "R1 size: ~${FILE_SIZE_MB} MB"

THREAD_ARGS=""
if [[ -n "$THREADS" ]]; then
  THREAD_ARGS="--threads $THREADS"
fi

echo ""
echo "=== Running hyperfine benchmarks ==="
echo "(sample_size=10, warmup=1)"
echo ""

hyperfine \
  --warmup 1 \
  --runs 10 \
  --export-markdown "$TMPDIR/results.md" \
  --export-json    "$TMPDIR/results.json" \
  --command-name "single-end (plain)" \
    "$BIN count --r1 $R1 --output $OUT $THREAD_ARGS" \
  --command-name "single-end (gzip)" \
    "$BIN count --r1 $R1_GZ --output $OUT $THREAD_ARGS" \
  --command-name "paired-end (plain)" \
    "$BIN count --r1 $R1 --r2 $R2 --output $OUT $THREAD_ARGS" \
  --command-name "single-end + trim-start" \
    "$BIN count --r1 $R1 --trim-start ACGT --output $OUT $THREAD_ARGS"

echo ""
echo "=== Results ==="
cat "$TMPDIR/results.md"

# Copy results to workspace
cp "$TMPDIR/results.md"   bench_results.md
cp "$TMPDIR/results.json" bench_results.json
echo ""
echo "Results saved to bench_results.md and bench_results.json"

# Thread scaling (only if not already fixed)
if [[ -z "$THREADS" ]]; then
  echo ""
  echo "=== Thread scaling (single-end, $READS reads) ==="
  THREAD_CMDS=()
  for T in 1 2 4 8; do
    THREAD_CMDS+=("--command-name" "threads=$T" \
                  "$BIN count --r1 $R1 --output $OUT --threads $T")
  done
  hyperfine \
    --warmup 1 \
    --runs 5 \
    --export-markdown bench_scaling.md \
    "${THREAD_CMDS[@]}"
  echo ""
  echo "Scaling results saved to bench_scaling.md"
fi
