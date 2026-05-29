use ahash::AHashMap;
use clap::{Parser, Subcommand};
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};

/// Number of FASTQ records per batch handed off to the worker pool.
const BATCH_SIZE: usize = 50_000;


////////////////////////////////////////////////////////////////////////////////
// CLI
////////////////////////////////////////////////////////////////////////////////

#[derive(Parser, Debug)]
#[command(author, version, about = "Count FASTQ read sequences by identity")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Count read sequences across one (single-end) or two (paired-end) FASTQs.
    Count(CountArgs),
    /// Match a counts TSV against a predefined-sequences TSV.
    Match(MatchArgs),
}

#[derive(Parser, Debug)]
struct CountArgs {
    /// R1 FASTQ file (plain or gzip).
    #[arg(short = '1', long)]
    r1: String,

    /// R2 FASTQ file (plain or gzip). Omit for single-end mode.
    #[arg(short = '2', long)]
    r2: Option<String>,

    /// Trim read from the first occurrence of this k-mer (inclusive).
    #[arg(long)]
    trim_start: Option<String>,

    /// Trim read up to (exclusive) the last occurrence of this k-mer.
    #[arg(long)]
    trim_stop: Option<String>,

    /// Keep at most this many bases after adapter trimming.
    #[arg(long)]
    trim_length: Option<usize>,

    /// Output TSV file. Writes to stdout when omitted.
    #[arg(short, long)]
    output: Option<String>,

    /// Number of worker threads for counting (default: all logical CPUs).
    #[arg(short = 't', long)]
    threads: Option<usize>,
}

#[derive(Parser, Debug)]
struct MatchArgs {
    /// Counts TSV produced by the `count` subcommand.
    #[arg(long)]
    counts: String,

    /// Predefined-sequences TSV.
    #[arg(long)]
    predefined: String,

    /// Column in the predefined TSV holding the R1 sequence. Default: "Sequence".
    #[arg(long, default_value = "Sequence")]
    seq_col: String,

    /// Column in the predefined TSV holding the R2 sequence (paired mode).
    #[arg(long)]
    r2_col: Option<String>,

    /// Column in the predefined TSV holding the sequence identifier. Default: "Name".
    #[arg(long, default_value = "Name")]
    id_col: String,

    /// Output TSV for matched sequences.
    #[arg(short, long)]
    output: Option<String>,

    /// Output TSV for unmatched (counted but not predefined) sequences.
    #[arg(long)]
    unmatched: Option<String>,
}


////////////////////////////////////////////////////////////////////////////////
// FASTQ I/O
////////////////////////////////////////////////////////////////////////////////

fn open_reader(path: &str) -> io::Result<BufReader<Box<dyn Read>>> {
    let file = File::open(path)?;
    let inner: Box<dyn Read> = if path.ends_with(".gz") {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Ok(BufReader::with_capacity(1 << 20, inner))
}

struct FastqIter<R: BufRead> {
    reader: R,
    line: String,
}

impl<R: BufRead> FastqIter<R> {
    fn new(reader: R) -> Self {
        FastqIter {
            reader,
            line: String::with_capacity(512),
        }
    }

    fn next_record(&mut self) -> Option<(String, String)> {
        // @header
        self.line.clear();
        match self.reader.read_line(&mut self.line) {
            Ok(0) => return None,
            Ok(_) => {}
            Err(e) => panic!("FASTQ read error: {e}"),
        }
        if !self.line.starts_with('@') {
            panic!("Expected '@' header, got: {}", self.line.trim());
        }
        let name = self.line[1..]
            .trim()
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_owned();

        // sequence
        self.line.clear();
        self.reader.read_line(&mut self.line).expect("Truncated FASTQ: missing sequence");
        let seq = self.line.trim_end().to_owned();

        // '+' separator
        self.line.clear();
        self.reader.read_line(&mut self.line).expect("Truncated FASTQ: missing '+'");

        // quality (discard)
        self.line.clear();
        self.reader.read_line(&mut self.line).expect("Truncated FASTQ: missing quality");

        Some((name, seq))
    }
}


////////////////////////////////////////////////////////////////////////////////
// Trimming
////////////////////////////////////////////////////////////////////////////////

/// Trim a sequence by optional start/stop k-mers and optional max length.
/// Returns empty string if trim_start is given but not found (invalid read).
fn trim_sequence<'a>(
    seq: &'a str,
    trim_start: Option<&str>,
    trim_stop: Option<&str>,
    length: Option<usize>,
) -> &'a str {
    let mut s = seq;

    if let Some(kmer) = trim_start {
        match s.find(kmer) {
            Some(pos) => s = &s[pos..],
            None => return "",
        }
    }

    if let Some(kmer) = trim_stop {
        if let Some(pos) = s.rfind(kmer) {
            s = &s[..pos];
        }
    }

    if let Some(n) = length {
        let end = s
            .char_indices()
            .nth(n)
            .map(|(i, _)| i)
            .unwrap_or(s.len());
        s = &s[..end];
    }

    s
}


////////////////////////////////////////////////////////////////////////////////
// Count data structure
////////////////////////////////////////////////////////////////////////////////

struct CountEntry {
    r1: String,
    r2: Option<String>,
    count: u64,
    r1_name: String,
    r2_name: Option<String>,
}


////////////////////////////////////////////////////////////////////////////////
// Paired counting
////////////////////////////////////////////////////////////////////////////////

fn count_paired(
    r1_path: &str,
    r2_path: &str,
    trim_start: Option<&str>,
    trim_stop: Option<&str>,
    trim_length: Option<usize>,
) -> io::Result<Vec<CountEntry>> {
    eprintln!("Counting paired reads: {r1_path} + {r2_path}");

    // --- Reader thread: reads R1+R2 in lockstep, sends batches via channel ---
    let (tx, rx) = std::sync::mpsc::sync_channel::<Vec<(String, String, String, String)>>(16);
    let r1_owned = r1_path.to_owned();
    let r2_owned = r2_path.to_owned();

    let reader = std::thread::spawn(move || -> io::Result<()> {
        let mut iter1 = FastqIter::new(open_reader(&r1_owned)?);
        let mut iter2 = FastqIter::new(open_reader(&r2_owned)?);
        let mut batch: Vec<(String, String, String, String)> =
            Vec::with_capacity(BATCH_SIZE);
        loop {
            match (iter1.next_record(), iter2.next_record()) {
                (Some((n1, s1)), Some((n2, s2))) => {
                    batch.push((n1, s1, n2, s2));
                    if batch.len() >= BATCH_SIZE {
                        if tx
                            .send(std::mem::replace(
                                &mut batch,
                                Vec::with_capacity(BATCH_SIZE),
                            ))
                            .is_err()
                        {
                            break;
                        }
                    }
                }
                (None, None) => break,
                _ => panic!("R1 and R2 have different numbers of records"),
            }
        }
        if !batch.is_empty() {
            let _ = tx.send(batch);
        }
        Ok(())
    });

    let batches: Vec<Vec<(String, String, String, String)>> = rx.iter().collect();
    reader.join().unwrap()?;

    // --- Worker pool: parallel fold over batches, then merge ---
    let trim_start_s = trim_start.map(str::to_owned);
    let trim_stop_s = trim_stop.map(str::to_owned);

    type PairedMap = AHashMap<(String, String), u64>;
    type PairedNames = AHashMap<(String, String), (String, String)>;

    let (counts, mut names): (PairedMap, PairedNames) = batches
        .into_par_iter()
        .fold(
            || (PairedMap::new(), PairedNames::new()),
            |(mut counts, mut names), batch| {
                for (n1, s1, n2, s2) in batch {
                    let t1 = trim_sequence(
                        &s1,
                        trim_start_s.as_deref(),
                        trim_stop_s.as_deref(),
                        trim_length,
                    );
                    let t2 = trim_sequence(
                        &s2,
                        trim_start_s.as_deref(),
                        trim_stop_s.as_deref(),
                        trim_length,
                    );
                    if t1.is_empty() || t2.is_empty() {
                        continue;
                    }
                    let key = (t1.to_owned(), t2.to_owned());
                    *counts.entry(key.clone()).or_insert(0) += 1;
                    names.entry(key).or_insert((n1, n2));
                }
                (counts, names)
            },
        )
        .reduce(
            || (PairedMap::new(), PairedNames::new()),
            |(mut ca, mut na), (cb, nb)| {
                for (k, v) in cb {
                    *ca.entry(k).or_insert(0) += v;
                }
                for (k, v) in nb {
                    na.entry(k).or_insert(v);
                }
                (ca, na)
            },
        );

    let mut result: Vec<CountEntry> = counts
        .into_iter()
        .map(|(key, count)| {
            let (n1, n2) = names.remove(&key).unwrap_or_default();
            CountEntry {
                r1: key.0,
                r2: Some(key.1),
                count,
                r1_name: n1,
                r2_name: Some(n2),
            }
        })
        .collect();

    result.sort_unstable_by(|a, b| b.count.cmp(&a.count));
    Ok(result)
}


////////////////////////////////////////////////////////////////////////////////
// Single-end counting
////////////////////////////////////////////////////////////////////////////////

fn count_single(
    r1_path: &str,
    trim_start: Option<&str>,
    trim_stop: Option<&str>,
    trim_length: Option<usize>,
) -> io::Result<Vec<CountEntry>> {
    eprintln!("Counting single-end reads: {r1_path}");

    // --- Reader thread: reads records and sends batches via channel ---
    let (tx, rx) = std::sync::mpsc::sync_channel::<Vec<(String, String)>>(16);
    let r1_owned = r1_path.to_owned();

    let reader = std::thread::spawn(move || -> io::Result<()> {
        let mut iter = FastqIter::new(open_reader(&r1_owned)?);
        let mut batch: Vec<(String, String)> = Vec::with_capacity(BATCH_SIZE);
        while let Some(rec) = iter.next_record() {
            batch.push(rec);
            if batch.len() >= BATCH_SIZE {
                if tx
                    .send(std::mem::replace(
                        &mut batch,
                        Vec::with_capacity(BATCH_SIZE),
                    ))
                    .is_err()
                {
                    break;
                }
            }
        }
        if !batch.is_empty() {
            let _ = tx.send(batch);
        }
        Ok(())
    });

    let batches: Vec<Vec<(String, String)>> = rx.iter().collect();
    reader.join().unwrap()?;

    // --- Worker pool: parallel fold over batches, then merge ---
    let trim_start_s = trim_start.map(str::to_owned);
    let trim_stop_s = trim_stop.map(str::to_owned);

    type SingleMap = AHashMap<String, u64>;
    type SingleNames = AHashMap<String, String>;

    let (counts, mut names): (SingleMap, SingleNames) = batches
        .into_par_iter()
        .fold(
            || (SingleMap::new(), SingleNames::new()),
            |(mut counts, mut names), batch| {
                for (name, seq) in batch {
                    let t = trim_sequence(
                        &seq,
                        trim_start_s.as_deref(),
                        trim_stop_s.as_deref(),
                        trim_length,
                    );
                    if t.is_empty() {
                        continue;
                    }
                    let key = t.to_owned();
                    *counts.entry(key.clone()).or_insert(0) += 1;
                    names.entry(key).or_insert(name);
                }
                (counts, names)
            },
        )
        .reduce(
            || (SingleMap::new(), SingleNames::new()),
            |(mut ca, mut na), (cb, nb)| {
                for (k, v) in cb {
                    *ca.entry(k).or_insert(0) += v;
                }
                for (k, v) in nb {
                    na.entry(k).or_insert(v);
                }
                (ca, na)
            },
        );

    let mut result: Vec<CountEntry> = counts
        .into_iter()
        .map(|(seq, count)| {
            let n = names.remove(&seq).unwrap_or_default();
            CountEntry {
                r1: seq,
                r2: None,
                count,
                r1_name: n,
                r2_name: None,
            }
        })
        .collect();

    result.sort_unstable_by(|a, b| b.count.cmp(&a.count));
    Ok(result)
}

////////////////////////////////////////////////////////////////////////////////
// TSV output
////////////////////////////////////////////////////////////////////////////////


fn write_count_tsv(mut w: Box<dyn Write>, entries: &[CountEntry]) -> io::Result<()> {
    let paired = entries.first().map_or(false, |e| e.r2.is_some());
    if paired {
        writeln!(w, "R1\tR2\tCount\tR1 Name\tR2 Name")?;
        for e in entries {
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}",
                e.r1,
                e.r2.as_deref().unwrap_or(""),
                e.count,
                e.r1_name,
                e.r2_name.as_deref().unwrap_or("")
            )?;
        }
    } else {
        writeln!(w, "R1\tCount\tR1 Name")?;
        for e in entries {
            writeln!(w, "{}\t{}\t{}", e.r1, e.count, e.r1_name)?;
        }
    }
    Ok(())
}

fn make_writer(path: &Option<String>) -> io::Result<Box<dyn Write>> {
    match path {
        Some(p) => Ok(Box::new(BufWriter::new(File::create(p)?))),
        None => Ok(Box::new(BufWriter::new(io::stdout().lock()))),
    }
}


////////////////////////////////////////////////////////////////////////////////
// Count subcommand runner
////////////////////////////////////////////////////////////////////////////////

fn run_count(args: &CountArgs) -> io::Result<()> {
    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .ok();
    }

    let trim_start = args.trim_start.as_deref();
    let trim_stop = args.trim_stop.as_deref();
    let trim_length = args.trim_length;

    let entries = if let Some(r2_path) = &args.r2 {
        count_paired(&args.r1, r2_path, trim_start, trim_stop, trim_length)?
    } else {
        count_single(&args.r1, trim_start, trim_stop, trim_length)?
    };

    let total: u64 = entries.iter().map(|e| e.count).sum();
    eprintln!(
        "Done. {} unique sequences, {} total reads.",
        entries.len(),
        total
    );

    write_count_tsv(make_writer(&args.output)?, &entries)?;

    if let Some(p) = &args.output {
        eprintln!("Written to {p}");
    }
    Ok(())
}


////////////////////////////////////////////////////////////////////////////////
// TSV reader
////////////////////////////////////////////////////////////////////////////////

fn read_tsv(path: &str) -> io::Result<(Vec<String>, Vec<Vec<String>>)> {
    let reader = BufReader::new(File::open(path)?);
    let mut lines = reader.lines();

    let header: Vec<String> = lines
        .next()
        .ok_or_else(|| io::Error::new(io::ErrorKind::UnexpectedEof, "Empty TSV"))??
        .split('\t')
        .map(str::to_owned)
        .collect();

    let rows: Vec<Vec<String>> = lines
        .map(|l| l.map(|s| s.split('\t').map(str::to_owned).collect()))
        .collect::<io::Result<_>>()?;

    Ok((header, rows))
}

fn col_idx(header: &[String], name: &str) -> io::Result<usize> {
    header.iter().position(|h| h == name).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Column '{name}' not found in TSV header"),
        )
    })
}


////////////////////////////////////////////////////////////////////////////////
// Match subcommand runner
////////////////////////////////////////////////////////////////////////////////

fn run_match(args: &MatchArgs) -> io::Result<()> {
    let (count_header, count_rows) = read_tsv(&args.counts)?;
    let (pred_header, pred_rows) = read_tsv(&args.predefined)?;

    let c_r1 = col_idx(&count_header, "R1")?;
    let c_r2 = count_header.iter().position(|h| h == "R2");
    let c_count = col_idx(&count_header, "Count")?;

    let p_seq = col_idx(&pred_header, &args.seq_col)?;
    let p_r2 = args
        .r2_col
        .as_deref()
        .map(|col| col_idx(&pred_header, col))
        .transpose()?;

    let total_reads: u64 = count_rows
        .iter()
        .filter_map(|r| r.get(c_count).and_then(|v| v.parse::<u64>().ok()))
        .sum();

    // Build lookup: key -> count
    let mut count_map: AHashMap<String, u64> = AHashMap::new();
    for row in &count_rows {
        let key = match (c_r2, p_r2) {
            (Some(ci), Some(_)) => format!(
                "{}\t{}",
                row.get(c_r1).map_or("", String::as_str),
                row.get(ci).map_or("", String::as_str)
            ),
            _ => row.get(c_r1).map_or("", String::as_str).to_owned(),
        };
        if let Some(v) = row.get(c_count).and_then(|v| v.parse::<u64>().ok()) {
            count_map.insert(key, v);
        }
    }

    // Write matched
    let mut matched_writer = make_writer(&args.output)?;
    {
        let mut out_header = pred_header.clone();
        out_header.extend(["Read Count".to_owned(), "Frequency".to_owned()]);
        writeln!(matched_writer, "{}", out_header.join("\t"))?;
    }

    let mut matched_keys: std::collections::HashSet<String> = std::collections::HashSet::new();

    for row in &pred_rows {
        let key = match p_r2 {
            Some(pi) => format!(
                "{}\t{}",
                row.get(p_seq).map_or("", String::as_str),
                row.get(pi).map_or("", String::as_str)
            ),
            None => row.get(p_seq).map_or("", String::as_str).to_owned(),
        };
        let count = count_map.get(&key).copied().unwrap_or(0);
        let freq = if total_reads > 0 {
            count as f64 / total_reads as f64
        } else {
            0.0
        };
        writeln!(
            matched_writer,
            "{}\t{}\t{:.6}",
            row.join("\t"),
            count,
            freq
        )?;
        matched_keys.insert(key);
    }

    // Write unmatched
    if let Some(unmatched_path) = &args.unmatched {
        let mut uw = make_writer(&Some(unmatched_path.clone()))?;
        writeln!(uw, "{}", count_header.join("\t"))?;
        for row in &count_rows {
            let key = match (c_r2, p_r2) {
                (Some(ci), Some(_)) => format!(
                    "{}\t{}",
                    row.get(c_r1).map_or("", String::as_str),
                    row.get(ci).map_or("", String::as_str)
                ),
                _ => row.get(c_r1).map_or("", String::as_str).to_owned(),
            };
            if !matched_keys.contains(&key) {
                writeln!(uw, "{}", row.join("\t"))?;
            }
        }
        eprintln!("Unmatched written to {unmatched_path}");
    }

    eprintln!(
        "Matched {} predefined sequences against {} counted entries.",
        pred_rows.len(),
        count_rows.len()
    );
    Ok(())
}


////////////////////////////////////////////////////////////////////////////////
// Entry point
////////////////////////////////////////////////////////////////////////////////

fn main() -> io::Result<()> {
    let cli = Cli::parse();
    match &cli.command {
        Command::Count(args) => run_count(args),
        Command::Match(args) => run_match(args),
    }
}