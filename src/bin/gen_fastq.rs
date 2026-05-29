/// Generate synthetic FASTQ files for benchmarking.
///
/// Usage:
///   gen_fastq --reads 1000000 --read-len 150 --num-seqs 500 --output r1.fastq
///   gen_fastq --reads 1000000 --read-len 150 --num-seqs 500 --output r1.fastq --output-r2 r2.fastq
///   gen_fastq --reads 1000000 --read-len 150 --num-seqs 500 --output r1.fastq.gz

/// cargo run --release --bin gen_fastq --reads 5000000 --read-len 150 --num-seqs 10000 --output r1.fastq --output-r2 r2.fastq

// for crierion: cargo bench
// # Reports: target/criterion/index.html

// with hyperfine:
// (cargo install hyperfine)
// bash scripts/bench.sh --reads 10000000 --threads 8
// results in bench_results.md und bench_scaling.md
use clap::Parser;
use flate2::write::GzEncoder;
use flate2::Compression;
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use std::fs::File;
use std::io::{BufWriter, Write};

#[derive(Parser, Debug)]
#[command(about = "Generate synthetic FASTQ files for benchmarking")]
struct Args {
    /// Number of reads to generate.
    #[arg(short = 'n', long, default_value = "1000000")]
    reads: usize,

    /// Read length (bases).
    #[arg(short = 'l', long, default_value = "150")]
    read_len: usize,

    /// Number of distinct sequences in the pool (controls uniqueness).
    #[arg(short = 's', long, default_value = "1000")]
    num_seqs: usize,

    /// Output R1 FASTQ file (use .fastq.gz for gzip output).
    #[arg(short = 'o', long)]
    output: String,

    /// Output R2 FASTQ file (optional, for paired-end output).
    #[arg(long)]
    output_r2: Option<String>,

    /// RNG seed for reproducibility.
    #[arg(long, default_value = "42")]
    seed: u64,
}

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
const QUAL_CHAR: u8 = b'I'; // Phred 40

fn make_writer(path: &str) -> Box<dyn Write> {
    let file = File::create(path).expect("Cannot create output file");
    if path.ends_with(".gz") {
        Box::new(BufWriter::new(GzEncoder::new(file, Compression::default())))
    } else {
        Box::new(BufWriter::new(file))
    }
}

fn gen_sequence(rng: &mut SmallRng, len: usize) -> Vec<u8> {
    (0..len).map(|_| BASES[rng.gen_range(0..4)]).collect()
}

fn main() {
    let args = Args::parse();

    let mut rng = SmallRng::seed_from_u64(args.seed);

    // Pre-generate the pool of distinct sequences.
    eprintln!(
        "Generating pool of {} distinct sequences (len={})…",
        args.num_seqs, args.read_len
    );
    let pool_r1: Vec<Vec<u8>> = (0..args.num_seqs)
        .map(|_| gen_sequence(&mut rng, args.read_len))
        .collect();
    let pool_r2: Vec<Vec<u8>> = (0..args.num_seqs)
        .map(|_| gen_sequence(&mut rng, args.read_len))
        .collect();
    let qual: Vec<u8> = vec![QUAL_CHAR; args.read_len];

    let paired = args.output_r2.is_some();

    eprintln!(
        "Writing {} reads to {} {}…",
        args.reads,
        args.output,
        if paired {
            format!("+ {}", args.output_r2.as_deref().unwrap())
        } else {
            String::new()
        }
    );

    let mut w1 = make_writer(&args.output);
    let mut w2 = args.output_r2.as_deref().map(make_writer);

    for i in 0..args.reads {
        let idx = rng.gen_range(0..args.num_seqs);

        // R1
        writeln!(w1, "@read{i}").unwrap();
        w1.write_all(&pool_r1[idx]).unwrap();
        writeln!(w1).unwrap();
        writeln!(w1, "+").unwrap();
        w1.write_all(&qual).unwrap();
        writeln!(w1).unwrap();

        // R2 (optional)
        if let Some(ref mut w) = w2 {
            writeln!(w, "@read{i}").unwrap();
            w.write_all(&pool_r2[idx]).unwrap();
            writeln!(w).unwrap();
            writeln!(w, "+").unwrap();
            w.write_all(&qual).unwrap();
            writeln!(w).unwrap();
        }
    }

    eprintln!("Done.");
}
