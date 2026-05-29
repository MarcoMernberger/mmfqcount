//! Criterion benchmarks for mmfqcount.
//!
//! Each benchmark generates synthetic FASTQ data in a temp directory (setup
//! phase, excluded from timing), then times the `mmfqcount count` binary
//! invocations using `iter_custom` so the full wall-clock time per call
//! (I/O + parsing + trimming + hashing) is captured.
//!
//! Run with:
//!   cargo bench
//!
//! HTML reports land in target/criterion/.

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use std::fs;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{Duration, Instant};

// Path to the compiled binary, injected by Cargo at build time.
const BIN: &str = env!("CARGO_BIN_EXE_mmfqcount");

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
const QUAL_CHAR: u8 = b'I';

/// Generate a plain-text FASTQ with `n_reads` reads of length `read_len`,
/// drawn from a pool of `num_seqs` distinct sequences.
fn write_fastq(path: &Path, n_reads: usize, read_len: usize, num_seqs: usize, seed: u64) {
    let mut rng = SmallRng::seed_from_u64(seed);
    let pool: Vec<Vec<u8>> = (0..num_seqs)
        .map(|_| (0..read_len).map(|_| BASES[rng.gen_range(0..4)]).collect())
        .collect();
    let qual: Vec<u8> = vec![QUAL_CHAR; read_len];

    let mut w = BufWriter::new(fs::File::create(path).unwrap());
    for i in 0..n_reads {
        let idx = rng.gen_range(0..num_seqs);
        writeln!(w, "@r{i}").unwrap();
        w.write_all(&pool[idx]).unwrap();
        writeln!(w).unwrap();
        writeln!(w, "+").unwrap();
        w.write_all(&qual).unwrap();
        writeln!(w).unwrap();
    }
}

/// Run the binary and return wall-clock duration.
fn run_count_single(r1: &Path, out: &Path) -> Duration {
    let start = Instant::now();
    let status = Command::new(BIN)
        .args([
            "count",
            "--r1",
            r1.to_str().unwrap(),
            "--output",
            out.to_str().unwrap(),
        ])
        .stderr(std::process::Stdio::null())
        .status()
        .expect("failed to run mmfqcount");
    let elapsed = start.elapsed();
    assert!(status.success(), "mmfqcount exited non-zero");
    elapsed
}

fn run_count_paired(r1: &Path, r2: &Path, out: &Path) -> Duration {
    let start = Instant::now();
    let status = Command::new(BIN)
        .args([
            "count",
            "--r1",
            r1.to_str().unwrap(),
            "--r2",
            r2.to_str().unwrap(),
            "--output",
            out.to_str().unwrap(),
        ])
        .stderr(std::process::Stdio::null())
        .status()
        .expect("failed to run mmfqcount");
    let elapsed = start.elapsed();
    assert!(status.success(), "mmfqcount exited non-zero");
    elapsed
}

fn run_count_single_trim(r1: &Path, out: &Path, trim_start: &str) -> Duration {
    let start = Instant::now();
    let status = Command::new(BIN)
        .args([
            "count",
            "--r1",
            r1.to_str().unwrap(),
            "--trim-start",
            trim_start,
            "--output",
            out.to_str().unwrap(),
        ])
        .stderr(std::process::Stdio::null())
        .status()
        .expect("failed to run mmfqcount");
    let elapsed = start.elapsed();
    assert!(status.success(), "mmfqcount exited non-zero");
    elapsed
}

// ---------------------------------------------------------------------------
// Benchmark: single-end, varying read counts
// ---------------------------------------------------------------------------
fn bench_single_end(c: &mut Criterion) {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.tsv");

    let mut group = c.benchmark_group("count_single");

    for &n_reads in &[100_000usize, 500_000, 1_000_000] {
        let r1: PathBuf = dir.path().join(format!("r1_{n_reads}.fastq"));
        write_fastq(&r1, n_reads, 150, 1000, 42);

        group.throughput(Throughput::Elements(n_reads as u64));
        group.bench_with_input(
            BenchmarkId::new("reads", n_reads),
            &n_reads,
            |b, _| {
                b.iter_custom(|iters| {
                    (0..iters).map(|_| run_count_single(&r1, &out)).sum()
                });
            },
        );
    }
    group.finish();
}

// ---------------------------------------------------------------------------
// Benchmark: paired-end, varying read counts
// ---------------------------------------------------------------------------
fn bench_paired_end(c: &mut Criterion) {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.tsv");

    let mut group = c.benchmark_group("count_paired");

    for &n_reads in &[100_000usize, 500_000, 1_000_000] {
        let r1 = dir.path().join(format!("r1_{n_reads}.fastq"));
        let r2 = dir.path().join(format!("r2_{n_reads}.fastq"));
        write_fastq(&r1, n_reads, 150, 1000, 42);
        write_fastq(&r2, n_reads, 150, 1000, 99);

        group.throughput(Throughput::Elements(n_reads as u64));
        group.bench_with_input(
            BenchmarkId::new("reads", n_reads),
            &n_reads,
            |b, _| {
                b.iter_custom(|iters| {
                    (0..iters)
                        .map(|_| run_count_paired(&r1, &r2, &out))
                        .sum()
                });
            },
        );
    }
    group.finish();
}

// ---------------------------------------------------------------------------
// Benchmark: single-end + trim-start, varying read lengths
// ---------------------------------------------------------------------------
fn bench_trim_start(c: &mut Criterion) {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.tsv");
    const N: usize = 500_000;
    const ADAPTER: &str = "ACGT";

    let mut group = c.benchmark_group("count_single_trim_start");

    for &read_len in &[50usize, 100, 150, 250] {
        // Build reads that always have the adapter at position 0.
        let r1 = dir.path().join(format!("r1_trim_{read_len}.fastq"));
        {
            let mut rng = SmallRng::seed_from_u64(7);
            let pool: Vec<Vec<u8>> = (0..500)
                .map(|_| {
                    let mut seq = ADAPTER.as_bytes().to_vec();
                    seq.extend((0..read_len - ADAPTER.len()).map(|_| BASES[rng.gen_range(0..4)]));
                    seq
                })
                .collect();
            let qual = vec![QUAL_CHAR; read_len];
            let mut w = BufWriter::new(fs::File::create(&r1).unwrap());
            for i in 0..N {
                let idx = rng.gen_range(0..pool.len());
                writeln!(w, "@r{i}").unwrap();
                w.write_all(&pool[idx]).unwrap();
                writeln!(w).unwrap();
                writeln!(w, "+").unwrap();
                w.write_all(&qual).unwrap();
                writeln!(w).unwrap();
            }
        }

        group.throughput(Throughput::Elements(N as u64));
        group.bench_with_input(
            BenchmarkId::new("read_len", read_len),
            &read_len,
            |b, _| {
                b.iter_custom(|iters| {
                    (0..iters)
                        .map(|_| run_count_single_trim(&r1, &out, ADAPTER))
                        .sum()
                });
            },
        );
    }
    group.finish();
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(10);
    targets = bench_single_end, bench_paired_end, bench_trim_start
}
criterion_main!(benches);
