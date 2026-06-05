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
    /// Print a machine-readable JSON description of all subcommands and their parameters.
    Params(ParamsArgs),
}

#[derive(Parser, Debug)]
struct ParamsArgs {
    /// Write JSON to this file instead of stdout.
    #[arg(short, long)]
    output: Option<String>,
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

    /// Split counts by this tag in read names (e.g. "sgRNAid").
    /// The tag must appear as |TAG=VALUE| in the read name.
    /// Produces one 'Read Count (TAG=VALUE)' + 'Frequency (TAG=VALUE)' column pair per unique value.
    #[arg(long)]
    split_by: Option<String>,
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
    /// Per-tag-value counts; non-empty only when --split-by is used.
    tag_counts: AHashMap<String, u64>,
    r1_name: String,
    r2_name: Option<String>,
}

////////////////////////////////////////////////////////////////////////////////
// Tag extraction
////////////////////////////////////////////////////////////////////////////////

/// Extract the value of `tag` from a read name formatted as
/// `....|TAG=VALUE|...` (pipe-separated key=value pairs).
/// Returns `""` when the tag is absent.
fn extract_tag<'a>(name: &'a str, tag: &str) -> &'a str {
    // Build needle "|TAG=" to find the value start.
    let needle_len = tag.len() + 2; // '|' + tag + '='
    let bytes = name.as_bytes();
    let tag_bytes = tag.as_bytes();
    let mut i = 0;
    while i + needle_len <= bytes.len() {
        if bytes[i] == b'|'
            && bytes[i + 1..i + 1 + tag_bytes.len()] == *tag_bytes
            && bytes[i + 1 + tag_bytes.len()] == b'='
        {
            let val_start = i + needle_len;
            let val_end = bytes[val_start..]
                .iter()
                .position(|&b| b == b'|')
                .map(|p| val_start + p)
                .unwrap_or(name.len());
            return &name[val_start..val_end];
        }
        i += 1;
    }
    ""
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
    split_by: Option<&str>,
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
    let trim_stop_s  = trim_stop.map(str::to_owned);
    let split_by_s   = split_by.map(str::to_owned);

    // key: (t1, t2) -> tag_val -> count  (tag_val = "" when no split_by)
    type TagMap   = AHashMap<(String, String), AHashMap<String, u64>>;
    type NamesMap = AHashMap<(String, String), (String, String)>;

    let (tag_map, mut names): (TagMap, NamesMap) = batches
        .into_par_iter()
        .fold(
            || (TagMap::new(), NamesMap::new()),
            |(mut tag_map, mut names), batch| {
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
                    let tag_val = split_by_s
                        .as_deref()
                        .map(|tag| extract_tag(&n1, tag).to_owned())
                        .unwrap_or_default();
                    let key = (t1.to_owned(), t2.to_owned());
                    *tag_map.entry(key.clone()).or_default().entry(tag_val).or_insert(0) += 1;
                    names.entry(key).or_insert((n1, n2));
                }
                (tag_map, names)
            },
        )
        .reduce(
            || (TagMap::new(), NamesMap::new()),
            |(mut ta, mut na), (tb, nb)| {
                for (key, tvals) in tb {
                    let e = ta.entry(key).or_default();
                    for (tv, cnt) in tvals {
                        *e.entry(tv).or_insert(0) += cnt;
                    }
                }
                for (k, v) in nb {
                    na.entry(k).or_insert(v);
                }
                (ta, na)
            },
        );

    let has_split = split_by.is_some();
    let mut result: Vec<CountEntry> = tag_map
        .into_iter()
        .map(|(key, tvals)| {
            let count: u64 = tvals.values().sum();
            let (n1, n2) = names.remove(&key).unwrap_or_default();
            let tag_counts = if has_split { tvals } else { AHashMap::new() };
            CountEntry {
                r1: key.0,
                r2: Some(key.1),
                count,
                tag_counts,
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
    split_by: Option<&str>,
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
    let trim_stop_s  = trim_stop.map(str::to_owned);
    let split_by_s   = split_by.map(str::to_owned);

    // key: seq -> tag_val -> count  (tag_val = "" when no split_by)
    type TagMap   = AHashMap<String, AHashMap<String, u64>>;
    type NamesMap = AHashMap<String, String>;

    let (tag_map, mut names): (TagMap, NamesMap) = batches
        .into_par_iter()
        .fold(
            || (TagMap::new(), NamesMap::new()),
            |(mut tag_map, mut names), batch| {
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
                    let tag_val = split_by_s
                        .as_deref()
                        .map(|tag| extract_tag(&name, tag).to_owned())
                        .unwrap_or_default();
                    let key = t.to_owned();
                    *tag_map.entry(key.clone()).or_default().entry(tag_val).or_insert(0) += 1;
                    names.entry(key).or_insert(name);
                }
                (tag_map, names)
            },
        )
        .reduce(
            || (TagMap::new(), NamesMap::new()),
            |(mut ta, mut na), (tb, nb)| {
                for (seq, tvals) in tb {
                    let e = ta.entry(seq).or_default();
                    for (tv, cnt) in tvals {
                        *e.entry(tv).or_insert(0) += cnt;
                    }
                }
                for (k, v) in nb {
                    na.entry(k).or_insert(v);
                }
                (ta, na)
            },
        );

    let has_split = split_by.is_some();
    let mut result: Vec<CountEntry> = tag_map
        .into_iter()
        .map(|(seq, tvals)| {
            let count: u64 = tvals.values().sum();
            let n = names.remove(&seq).unwrap_or_default();
            let tag_counts = if has_split { tvals } else { AHashMap::new() };
            CountEntry {
                r1: seq,
                r2: None,
                count,
                tag_counts,
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


fn write_count_tsv(
    mut w: Box<dyn Write>,
    entries: &[CountEntry],
    split_by: Option<&str>,
    paired: bool,
) -> io::Result<()> {

    if let Some(tag_name) = split_by {
        // Collect all unique tag values across all entries, sorted.
        let mut all_tags: Vec<String> = {
            let mut set = std::collections::HashSet::new();
            for e in entries {
                for k in e.tag_counts.keys() {
                    set.insert(k.clone());
                }
            }
            let mut v: Vec<String> = set.into_iter().collect();
            v.sort();
            v
        };

        // Per-tag-value totals for frequency denominators.
        let mut tag_totals: AHashMap<&str, u64> = AHashMap::new();
        for e in entries {
            for (tv, &cnt) in &e.tag_counts {
                *tag_totals.entry(tv.as_str()).or_insert(0) += cnt;
            }
        }

        // Header
        let mut header: Vec<String> = if paired {
            vec!["R1".into(), "R2".into()]
        } else {
            vec!["R1".into()]
        };
        for tv in &all_tags {
            header.push(format!("Read Count ({}={})", tag_name, tv));
            header.push(format!("Frequency ({}={})", tag_name, tv));
        }
        if paired {
            header.extend(["R1 Name".into(), "R2 Name".into()]);
        } else {
            header.push("R1 Name".into());
        }
        writeln!(w, "{}", header.join("\t"))?;

        // Rows
        for e in entries {
            let mut row: Vec<String> = if paired {
                vec![e.r1.clone(), e.r2.as_deref().unwrap_or("").to_owned()]
            } else {
                vec![e.r1.clone()]
            };
            for tv in &all_tags {
                let cnt   = e.tag_counts.get(tv.as_str()).copied().unwrap_or(0);
                let total = tag_totals.get(tv.as_str()).copied().unwrap_or(0);
                let freq  = if total > 0 { cnt as f64 / total as f64 } else { 0.0 };
                row.push(cnt.to_string());
                row.push(format!("{:.6}", freq));
            }
            if paired {
                row.push(e.r1_name.clone());
                row.push(e.r2_name.as_deref().unwrap_or("").to_owned());
            } else {
                row.push(e.r1_name.clone());
            }
            writeln!(w, "{}", row.join("\t"))?;
        }
    } else {
        // No split_by — one Count + Frequency column for all reads.
        let total: u64 = entries.iter().map(|e| e.count).sum();
        if paired {
            writeln!(w, "R1\tR2\tCount\tFrequency\tR1 Name\tR2 Name")?;
            for e in entries {
                let freq = if total > 0 { e.count as f64 / total as f64 } else { 0.0 };
                writeln!(
                    w,
                    "{}\t{}\t{}\t{:.6}\t{}\t{}",
                    e.r1,
                    e.r2.as_deref().unwrap_or(""),
                    e.count,
                    freq,
                    e.r1_name,
                    e.r2_name.as_deref().unwrap_or("")
                )?;
            }
        } else {
            writeln!(w, "R1\tCount\tFrequency\tR1 Name")?;
            for e in entries {
                let freq = if total > 0 { e.count as f64 / total as f64 } else { 0.0 };
                writeln!(w, "{}\t{}\t{:.6}\t{}", e.r1, e.count, freq, e.r1_name)?;
            }
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
    let trim_stop  = args.trim_stop.as_deref();
    let trim_length = args.trim_length;
    let split_by   = args.split_by.as_deref();

    let entries = if let Some(r2_path) = &args.r2 {
        count_paired(&args.r1, r2_path, trim_start, trim_stop, trim_length, split_by)?
    } else {
        count_single(&args.r1, trim_start, trim_stop, trim_length, split_by)?
    };

    let total: u64 = entries.iter().map(|e| e.count).sum();
    eprintln!(
        "Done. {} unique sequences, {} total reads.",
        entries.len(),
        total
    );

    write_count_tsv(make_writer(&args.output)?, &entries, split_by, args.r2.is_some())?;

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
// Params subcommand
////////////////////////////////////////////////////////////////////////////////

fn run_params(args: &ParamsArgs) -> io::Result<()> {
    // Build the JSON description entirely with serde_json::json! macro —
    // no extra derive macros needed on the main structs.
    let doc = serde_json::json!({
        "tool": "mmfqcount",
        "version": env!("CARGO_PKG_VERSION"),
        "description": "Count FASTQ read sequences by identity, with optional adapter trimming and tag-based splitting.",
        "subcommands": [
            {
                "name": "count",
                "description": "Count read sequences across one (single-end) or two (paired-end) FASTQs.",
                "parameters": [
                    {
                        "name": "r1",
                        "flags": ["--r1", "-1"],
                        "type": "path",
                        "required": true,
                        "description": "R1 FASTQ file (plain or .gz)."
                    },
                    {
                        "name": "r2",
                        "flags": ["--r2", "-2"],
                        "type": "path",
                        "required": false,
                        "description": "R2 FASTQ file (plain or .gz). Omit for single-end mode."
                    },
                    {
                        "name": "output",
                        "flags": ["--output", "-o"],
                        "type": "path",
                        "required": false,
                        "description": "Output TSV file. Writes to stdout when omitted."
                    },
                    {
                        "name": "trim_start",
                        "flags": ["--trim-start"],
                        "type": "string",
                        "required": false,
                        "description": "Trim read from the first occurrence of this k-mer (inclusive). Reads lacking the k-mer are discarded."
                    },
                    {
                        "name": "trim_stop",
                        "flags": ["--trim-stop"],
                        "type": "string",
                        "required": false,
                        "description": "Trim read up to (exclusive) the last occurrence of this k-mer."
                    },
                    {
                        "name": "trim_length",
                        "flags": ["--trim-length"],
                        "type": "integer",
                        "required": false,
                        "description": "Keep at most this many bases after adapter trimming."
                    },
                    {
                        "name": "split_by",
                        "flags": ["--split-by"],
                        "type": "string",
                        "required": false,
                        "description": "Split counts by this tag in read names (e.g. 'sgRNAid'). The tag must appear as |TAG=VALUE| in the read name. Produces one 'Read Count (TAG=VALUE)' + 'Frequency (TAG=VALUE)' column pair per unique value."
                    },
                    {
                        "name": "threads",
                        "flags": ["--threads", "-t"],
                        "type": "integer",
                        "required": false,
                        "description": "Number of worker threads (default: all logical CPUs)."
                    }
                ],
                "output_format": {
                    "type": "tsv",
                    "modes": [
                        {
                            "condition": "single-end, no split_by",
                            "columns": ["R1", "Count", "Frequency", "R1 Name"]
                        },
                        {
                            "condition": "paired-end, no split_by",
                            "columns": ["R1", "R2", "Count", "Frequency", "R1 Name", "R2 Name"]
                        },
                        {
                            "condition": "single-end, with split_by TAG",
                            "columns": ["R1", "Read Count (TAG=<value>)", "Frequency (TAG=<value>)", "...", "R1 Name"],
                            "note": "One Read Count + Frequency column pair per unique tag value, sorted alphabetically."
                        },
                        {
                            "condition": "paired-end, with split_by TAG",
                            "columns": ["R1", "R2", "Read Count (TAG=<value>)", "Frequency (TAG=<value>)", "...", "R1 Name", "R2 Name"],
                            "note": "One Read Count + Frequency column pair per unique tag value, sorted alphabetically."
                        }
                    ]
                }
            },
            {
                "name": "match",
                "description": "Match a counts TSV (from 'count') against a predefined-sequences TSV.",
                "parameters": [
                    {
                        "name": "counts",
                        "flags": ["--counts"],
                        "type": "path",
                        "required": true,
                        "description": "Counts TSV produced by the 'count' subcommand."
                    },
                    {
                        "name": "predefined",
                        "flags": ["--predefined"],
                        "type": "path",
                        "required": true,
                        "description": "Predefined-sequences TSV."
                    },
                    {
                        "name": "seq_col",
                        "flags": ["--seq-col"],
                        "type": "string",
                        "required": false,
                        "default": "Sequence",
                        "description": "Column in the predefined TSV holding the R1 sequence."
                    },
                    {
                        "name": "r2_col",
                        "flags": ["--r2-col"],
                        "type": "string",
                        "required": false,
                        "description": "Column in the predefined TSV holding the R2 sequence (paired mode)."
                    },
                    {
                        "name": "id_col",
                        "flags": ["--id-col"],
                        "type": "string",
                        "required": false,
                        "default": "Name",
                        "description": "Column in the predefined TSV holding the sequence identifier."
                    },
                    {
                        "name": "output",
                        "flags": ["--output", "-o"],
                        "type": "path",
                        "required": false,
                        "description": "Output TSV for matched sequences. Writes to stdout when omitted."
                    },
                    {
                        "name": "unmatched",
                        "flags": ["--unmatched"],
                        "type": "path",
                        "required": false,
                        "description": "Output TSV for sequences that were counted but not found in the predefined list."
                    }
                ],
                "output_format": {
                    "type": "tsv",
                    "note": "All columns from the predefined TSV are preserved; 'Read Count' and 'Frequency' columns are appended.",
                    "appended_columns": ["Read Count", "Frequency"]
                }
            },
            {
                "name": "params",
                "description": "Print this machine-readable JSON parameter specification.",
                "parameters": [
                    {
                        "name": "output",
                        "flags": ["--output", "-o"],
                        "type": "path",
                        "required": false,
                        "description": "Write JSON to this file instead of stdout."
                    }
                ]
            }
        ]
    });

    let json_str = serde_json::to_string_pretty(&doc)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

    match &args.output {
        Some(path) => {
            let mut f = File::create(path)?;
            writeln!(f, "{json_str}")?;
            eprintln!("Params written to {path}");
        }
        None => println!("{json_str}"),
    }
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
        Command::Params(args) => run_params(args),
    }
}