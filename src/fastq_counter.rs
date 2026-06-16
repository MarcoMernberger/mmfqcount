use ahash::AHashMap;
use clap::{Parser, Subcommand};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::sync::Arc;

/// Explosion guard: panic if unique key count exceeds this.
const MAX_UNIQUE_KEYS: usize = 9_999_999;

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

    /// Split counts by this tag in read names (e.g. "sgRNAid").
    /// The tag must appear as |TAG=VALUE| in the read name.
    /// Annotation column will contain the tag value, or "UNKNOWN" if absent.
    #[arg(long)]
    split_by: Option<String>,

    /// Sort output by: count-desc (default), count-asc, sequence, annotation, none.
    #[arg(long, default_value = "count-desc")]
    sort_by: String,
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

/// How to sort the output TSV rows.
#[derive(Debug, Clone, Copy)]
enum SortBy {
    CountDesc,
    CountAsc,
    Sequence,
    Annotation,
    None,
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
        self.reader
            .read_line(&mut self.line)
            .expect("Truncated FASTQ: missing sequence");
        let seq = self.line.trim_end().to_owned();

        // '+' separator
        self.line.clear();
        self.reader
            .read_line(&mut self.line)
            .expect("Truncated FASTQ: missing '+'");

        // quality (discard)
        self.line.clear();
        self.reader
            .read_line(&mut self.line)
            .expect("Truncated FASTQ: missing quality");

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
        let end = s.char_indices().nth(n).map(|(i, _)| i).unwrap_or(s.len());
        s = &s[..end];
    }

    s
}

////////////////////////////////////////////////////////////////////////////////
// Count data structures
////////////////////////////////////////////////////////////////////////////////

/// Key: (R1_seq, R2_seq_or_empty, annotation).
/// Arc<str> avoids repeated heap allocation for frequently-seen sequences.
type Key = (Arc<str>, Arc<str>, Arc<str>);

struct EntryValue {
    count: u64,
    r1_name: Arc<str>,
    r2_name: Arc<str>,
}

type CountMap = AHashMap<Key, EntryValue>;

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

/// Derive the annotation string for a read.
/// - With `--split-by TAG`: returns the tag value, or `"UNKNOWN"` if absent.
/// - Without `--split-by`: returns `"ALL"`.
fn extract_annotation<'a>(name: &'a str, split_by: Option<&str>) -> &'a str {
    match split_by {
        Some(tag) => {
            let val = extract_tag(name, tag);
            if val.is_empty() {
                "UNKNOWN"
            } else {
                val
            }
        }
        None => "ALL",
    }
}

/// Build a `Key` from trimmed sequences and annotation.
fn build_key(t1: &str, t2: &str, ann: &str) -> Key {
    (Arc::from(t1), Arc::from(t2), Arc::from(ann))
}

/// Insert or increment an entry in `map`.
/// Panics if the number of unique keys would exceed `MAX_UNIQUE_KEYS`.
fn update_map(map: &mut CountMap, key: Key, r1_name: &str, r2_name: &str) {
    if map.len() >= MAX_UNIQUE_KEYS && !map.contains_key(&key) {
        panic!(
            "Key explosion: >{MAX_UNIQUE_KEYS} unique (R1, R2, Annotation) combinations detected"
        );
    }
    let entry = map.entry(key).or_insert_with(|| EntryValue {
        count: 0,
        r1_name: Arc::from(r1_name),
        r2_name: Arc::from(r2_name),
    });
    entry.count += 1;
}

/// Sanity-check a raw sequence from the FASTQ parser.
fn validate_record(seq: &str, name: &str, which: &str) -> bool {
    // is the record empty or freakishly long?
    if seq.is_empty() {
        eprintln!("Warning: Empty {which} sequence for read: {name}");
        return false;
    }
    if seq.len() >= 10_000 {
        eprintln!(
            "Warning: Suspiciously long {which} sequence ({} bp) for read: {name}",
            seq.len()
        );
        return false;
    }
    true
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
) -> io::Result<CountMap> {
    eprintln!("Counting paired reads: {r1_path} + {r2_path}");

    let mut map: CountMap = AHashMap::with_capacity(4096);
    let mut iter1 = FastqIter::new(open_reader(r1_path)?);
    let mut iter2 = FastqIter::new(open_reader(r2_path)?);

    loop {
        match (iter1.next_record(), iter2.next_record()) {
            (Some((n1, s1)), Some((n2, s2))) => {
                if !validate_record(&s1, &n1, "R1") || !validate_record(&s2, &n2, "R2") {
                    continue;
                }
                let t1 = trim_sequence(&s1, trim_start, trim_stop, trim_length);
                let t2 = trim_sequence(&s2, trim_start, trim_stop, trim_length);
                if t1.is_empty() || t2.is_empty() {
                    continue;
                }
                let ann = extract_annotation(&n1, split_by);
                let key = build_key(t1, t2, ann);
                update_map(&mut map, key, &n1, &n2);
            }
            (None, None) => break,
            _ => panic!("R1 and R2 have different numbers of records"),
        }
    }

    Ok(map)
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
) -> io::Result<CountMap> {
    eprintln!("Counting single-end reads: {r1_path}");

    let mut map: CountMap = AHashMap::with_capacity(4096);
    let mut iter = FastqIter::new(open_reader(r1_path)?);

    while let Some((name, seq)) = iter.next_record() {
        if !validate_record(&seq, &name, "R1") {
            continue;
        }
        let t1 = trim_sequence(&seq, trim_start, trim_stop, trim_length);
        if t1.is_empty() {
            continue;
        }
        let ann = extract_annotation(&name, split_by);
        let key = build_key(t1, "", ann);
        update_map(&mut map, key, &name, "");
    }

    Ok(map)
}

////////////////////////////////////////////////////////////////////////////////
// TSV output
////////////////////////////////////////////////////////////////////////////////

fn write_long_tsv(
    mut w: Box<dyn Write>,
    map: CountMap,
    paired: bool,
    sort_by: SortBy,
) -> io::Result<()> {
    let total: u64 = map.values().map(|e| e.count).sum();

    // Header
    if paired {
        writeln!(w, "R1\tR2\tAnnotation\tCount\tFrequency\tR1 Name\tR2 Name")?;
    } else {
        writeln!(w, "R1\tAnnotation\tCount\tFrequency\tR1 Name")?;
    }

    // Collect, filter, sort
    let mut entries: Vec<(Key, EntryValue)> =
        map.into_iter().filter(|(_, v)| v.count > 0).collect();

    match sort_by {
        SortBy::CountDesc => entries.sort_unstable_by(|a, b| b.1.count.cmp(&a.1.count)),
        SortBy::CountAsc => entries.sort_unstable_by(|a, b| a.1.count.cmp(&b.1.count)),
        SortBy::Sequence => entries.sort_unstable_by(|a, b| a.0 .0.cmp(&b.0 .0)),
        SortBy::Annotation => entries.sort_unstable_by(|a, b| a.0 .2.cmp(&b.0 .2)),
        SortBy::None => {}
    }

    // Rows — one per key
    for ((r1, r2, ann), val) in &entries {
        let freq = if total > 0 {
            val.count as f64 / total as f64
        } else {
            0.0
        };
        if paired {
            writeln!(
                w,
                "{r1}\t{r2}\t{ann}\t{}\t{freq:.6}\t{}\t{}",
                val.count, val.r1_name, val.r2_name
            )?;
        } else {
            writeln!(w, "{r1}\t{ann}\t{}\t{freq:.6}\t{}", val.count, val.r1_name)?;
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
    let trim_start = args.trim_start.as_deref();
    let trim_stop = args.trim_stop.as_deref();
    let trim_length = args.trim_length;
    let split_by = args.split_by.as_deref();

    let sort_by = match args.sort_by.as_str() {
        "count-asc" => SortBy::CountAsc,
        "sequence" => SortBy::Sequence,
        "annotation" => SortBy::Annotation,
        "none" => SortBy::None,
        _ => SortBy::CountDesc,
    };

    let (map, paired) = if let Some(r2_path) = &args.r2 {
        (
            count_paired(
                &args.r1,
                r2_path,
                trim_start,
                trim_stop,
                trim_length,
                split_by,
            )?,
            true,
        )
    } else {
        (
            count_single(&args.r1, trim_start, trim_stop, trim_length, split_by)?,
            false,
        )
    };

    let total: u64 = map.values().map(|e| e.count).sum();
    eprintln!(
        "Done. {} unique sequences, {} total reads.",
        map.len(),
        total
    );

    write_long_tsv(make_writer(&args.output)?, map, paired, sort_by)?;

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
        out_header.extend(["Count".to_owned(), "Frequency".to_owned()]);
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
        writeln!(matched_writer, "{}\t{}\t{:.6}", row.join("\t"), count, freq)?;
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
                        "description": "Split counts by this tag in read names (e.g. 'sgRNAid'). The tag must appear as |TAG=VALUE| in the read name. The Annotation column will contain the tag value, or 'UNKNOWN' if absent."
                    },
                    {
                        "name": "sort_by",
                        "flags": ["--sort-by"],
                        "type": "string",
                        "required": false,
                        "default": "count-desc",
                        "allowed_values": ["count-desc", "count-asc", "sequence", "annotation", "none"],
                        "description": "Sort output rows by: count-desc (default), count-asc, sequence (R1 alphabetically), annotation, or none (hash order)."
                    }
                ],
                "output_format": {
                    "type": "tsv",
                    "layout": "long",
                    "columns_single_end": ["R1", "Annotation", "Count", "Frequency", "R1 Name"],
                    "columns_paired_end": ["R1", "R2", "Annotation", "Count", "Frequency", "R1 Name", "R2 Name"],
                    "note": "One row per (R1, R2, Annotation) combination. Annotation is the split_by tag value, or 'ALL' if --split-by is not used. Only rows with Count > 0 are written."
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
                    "note": "All columns from the predefined TSV are preserved; 'Count' and 'Frequency' columns are appended.",
                    "appended_columns": ["Count", "Frequency"]
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

    let json_str =
        serde_json::to_string_pretty(&doc).map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

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
