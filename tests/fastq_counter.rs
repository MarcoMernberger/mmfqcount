//| Integration tests for fastq_counter.
//!
//! Each test writes temporary FASTQ / TSV files, runs the compiled binary via
//! `std::process::Command` and asserts on the output.

use std::fs;
use std::path::Path;

use tempfile;

use assert_cmd::Command;


////////////////////////////////////////////////////////////////////////////////
// Helpers
////////////////////////////////////////////////////////////////////////////////

/// Write a string to a file inside `dir` and return the full path as String.
fn write(dir: &Path, name: &str, content: &str) -> String {
    let p = dir.join(name);
    fs::write(&p, content).unwrap();
    p.to_string_lossy().into_owned()
}

/// Parse a TSV into (header: Vec<String>, rows: Vec<Vec<String>>).
fn parse_tsv(content: &str) -> (Vec<String>, Vec<Vec<String>>) {
    let mut lines = content.lines();
    let header: Vec<String> = lines
        .next()
        .unwrap()
        .split('\t')
        .map(str::to_owned)
        .collect();
    let rows: Vec<Vec<String>> = lines
        .filter(|l| !l.is_empty())
        .map(|l| l.split('\t').map(str::to_owned).collect())
        .collect();
    (header, rows)
}

/// Return the value of column `col` in `row` given `header`.
fn get<'a>(header: &[String], row: &'a [String], col: &str) -> &'a str {
    let idx = header.iter().position(|h| h == col).unwrap();
    row[idx].as_str()
}

////////////////////////////////////////////////////////////////////////////////
// FASTQ fixtures 
////////////////////////////////////////////////////////////////////////////////
// Format: @name\nSEQ\n+\nQUAL\n  (four lines per record)

/// Single-end FASTQ: 5 reads — AAAA×3, CCCC×1, GGGG×1
const SE_FASTQ: &str = "\
@read1\nAAAA\n+\nIIII\n\
@read2\nAAAA\n+\nIIII\n\
@read3\nCCCC\n+\nIIII\n\
@read4\nAAAA\n+\nIIII\n\
@read5\nGGGG\n+\nIIII\n";

/// R1 for paired-end: AAAA×3, CCCC×1, TTTT×1
const PE_R1: &str = "\
@r1\nAAAA\n+\nIIII\n\
@r2\nAAAA\n+\nIIII\n\
@r3\nCCCC\n+\nIIII\n\
@r4\nAAAA\n+\nIIII\n\
@r5\nTTTT\n+\nIIII\n";

/// R2 for paired-end: TTTT×3 paired with AAAA, GGGG×1, CCCC×1
const PE_R2: &str = "\
@r1\nTTTT\n+\nIIII\n\
@r2\nTTTT\n+\nIIII\n\
@r3\nGGGG\n+\nIIII\n\
@r4\nTTTT\n+\nIIII\n\
@r5\nCCCC\n+\nIIII\n";

/// Single-end FASTQ with adapter prefix: reads start with ADAPTER then payload.
/// ADAPTER = "ACGT"
/// Payloads: AAAA×3, CCCC×1, one read with no adapter (should be skipped).
const SE_ADAPTER_FASTQ: &str = "\
@a1\nACGTAAAA\n+\nIIIIIIII\n\
@a2\nACGTAAAA\n+\nIIIIIIII\n\
@a3\nACGTCCCC\n+\nIIIIIIII\n\
@a4\nACGTAAAA\n+\nIIIIIIII\n\
@a5\nNNNNAAAA\n+\nIIIIIIII\n";

////////////////////////////////////////////////////////////////////////////////
// Tests
////////////////////////////////////////////////////////////////////////////////


// Test 1: count — single-end, no trimming

#[test]
fn count_single_basic() {
    let dir = tempfile::tempdir().unwrap();
    let r1 = write(dir.path(), "r1.fastq", SE_FASTQ);
    let out = dir.path().join("out.tsv").to_string_lossy().into_owned();

Command::cargo_bin("mmfqcount").unwrap()
        .args(["count", "--r1", &r1, "--output", &out])
        .assert()
        .success();

    let content = fs::read_to_string(&out).unwrap();
    let (header, rows) = parse_tsv(&content);

    // Header check
    assert_eq!(header, ["R1", "Count", "Frequency", "R1 Name"]);

    // Correct number of unique sequences
    assert_eq!(rows.len(), 3, "expected 3 unique sequences");

    // Sorted descending by count
    let counts: Vec<u64> = rows
        .iter()
        .map(|r| get(&header, r, "Count").parse().unwrap())
        .collect();
    assert!(counts.windows(2).all(|w| w[0] >= w[1]), "not sorted descending");

    // AAAA must have count 3 and be first
    assert_eq!(get(&header, &rows[0], "R1"), "AAAA");
    assert_eq!(get(&header, &rows[0], "Count"), "3");

    // Example name must be a valid read name from the FASTQ
    let name = get(&header, &rows[0], "R1 Name");
    assert!(
        ["read1", "read2", "read4"].contains(&name),
        "unexpected read name: {name}"
    );
}

// Test 2: count — paired-end, no trimming 

#[test]
fn count_paired_basic() {
    let dir = tempfile::tempdir().unwrap();
    let r1 = write(dir.path(), "r1.fastq", PE_R1);
    let r2 = write(dir.path(), "r2.fastq", PE_R2);
    let out = dir.path().join("out.tsv").to_string_lossy().into_owned();

Command::cargo_bin("mmfqcount").unwrap()
        .args(["count", "--r1", &r1, "--r2", &r2, "--output", &out])
        .assert()
        .success();

    let content = fs::read_to_string(&out).unwrap();
    let (header, rows) = parse_tsv(&content);

    assert_eq!(header, ["R1", "R2", "Count", "Frequency", "R1 Name", "R2 Name"]);
    assert_eq!(rows.len(), 3, "expected 3 unique pairs");

    // Sorted descending
    let counts: Vec<u64> = rows
        .iter()
        .map(|r| get(&header, r, "Count").parse().unwrap())
        .collect();
    assert!(counts.windows(2).all(|w| w[0] >= w[1]));

    // Most common pair: (AAAA, TTTT) × 3
    assert_eq!(get(&header, &rows[0], "R1"), "AAAA");
    assert_eq!(get(&header, &rows[0], "R2"), "TTTT");
    assert_eq!(get(&header, &rows[0], "Count"), "3");
}

// Test 3: count — single-end, trim-start adapter 

#[test]
fn count_single_trim_start() {
    let dir = tempfile::tempdir().unwrap();
    let r1 = write(dir.path(), "r1.fastq", SE_ADAPTER_FASTQ);
    let out = dir.path().join("out.tsv").to_string_lossy().into_owned();

Command::cargo_bin("mmfqcount").unwrap()
        .args([
            "count",
            "--r1", &r1,
            "--trim-start", "ACGT",
            "--output", &out,
        ])
        .assert()
        .success();

    let content = fs::read_to_string(&out).unwrap();
    let (header, rows) = parse_tsv(&content);

    // Read 5 ("NNNNAAAA") has no "ACGT" adapter → must be skipped
    // Remaining: ACGTAAAA×3 → trimmed ACGTAAAA, ACGTCCCC×1 → ACGTCCCC
    // (trim_start keeps from the first occurrence of the kmer inclusive)
    assert_eq!(rows.len(), 2, "read without adapter must be skipped");

    // The adapter is included in the trimmed sequence
    assert_eq!(get(&header, &rows[0], "R1"), "ACGTAAAA");
    assert_eq!(get(&header, &rows[0], "Count"), "3");
}

// Test 4: count — single-end, trim-start + trim-length

#[test]
fn count_single_trim_start_and_length() {
    let dir = tempfile::tempdir().unwrap();
    let r1 = write(dir.path(), "r1.fastq", SE_ADAPTER_FASTQ);
    let out = dir.path().join("out.tsv").to_string_lossy().into_owned();

    // trim adapter, then keep 4 bases → strips the "ACGT" prefix leaving payload
Command::cargo_bin("mmfqcount").unwrap()
        .args([
            "count",
            "--r1", &r1,
            "--trim-start", "ACGT",
            "--trim-length", "8",   // keep full 8-base trimmed sequence
            "--output", &out,
        ])
        .assert()
        .success();

    let content = fs::read_to_string(&out).unwrap();
    let (header, rows) = parse_tsv(&content);

    // Still 2 unique trimmed sequences (adapter-less reads skipped)
    assert_eq!(rows.len(), 2);
    // trim-length=8 keeps full "ACGTAAAA" (8 chars)
    assert_eq!(get(&header, &rows[0], "R1"), "ACGTAAAA");
}

// Test 5: count — single-end, trim-stop

#[test]
fn count_single_trim_stop() {
    // Reads: AAAASTOP, CCCCSTOP — trim_stop="STOP" → AAAA, CCCC
    let fastq = "\
@s1\nAAAASTOP\n+\nIIIIIIII\n\
@s2\nAAAASTOP\n+\nIIIIIIII\n\
@s3\nCCCCSTOP\n+\nIIIIIIII\n";

    let dir = tempfile::tempdir().unwrap();
    let r1 = write(dir.path(), "r1.fastq", fastq);
    let out = dir.path().join("out.tsv").to_string_lossy().into_owned();

Command::cargo_bin("mmfqcount").unwrap()
        .args([
            "count",
            "--r1", &r1,
            "--trim-stop", "STOP",
            "--output", &out,
        ])
        .assert()
        .success();

    let content = fs::read_to_string(&out).unwrap();
    let (header, rows) = parse_tsv(&content);

    assert_eq!(rows.len(), 2);
    assert_eq!(get(&header, &rows[0], "R1"), "AAAA");
    assert_eq!(get(&header, &rows[0], "Count"), "2");
    assert_eq!(get(&header, &rows[1], "R1"), "CCCC");
}

// Test 6: match — single-end

#[test]
fn match_single_basic() {
    // counts.tsv produced by a previous `count` run
    let counts_tsv = "R1\tCount\tR1 Name\nAAAA\t3\tread1\nCCCC\t1\tread3\nGGGG\t1\tread5\n";
    // predefined: only AAAA and GGGG
    let pred_tsv = "Name\tSequence\nseq_A\tAAAA\nseq_G\tGGGG\n";

    let dir = tempfile::tempdir().unwrap();
    let counts = write(dir.path(), "counts.tsv", counts_tsv);
    let predefined = write(dir.path(), "pred.tsv", pred_tsv);
    let matched = dir.path().join("matched.tsv").to_string_lossy().into_owned();
    let unmatched = dir.path().join("unmatched.tsv").to_string_lossy().into_owned();

Command::cargo_bin("mmfqcount").unwrap()
        .args([
            "match",
            "--counts",     &counts,
            "--predefined", &predefined,
            "--seq-col",    "Sequence",
            "--id-col",     "Name",
            "--output",     &matched,
            "--unmatched",  &unmatched,
        ])
        .assert()
        .success();

    // Matched output
    let mc = fs::read_to_string(&matched).unwrap();
    let (mh, mr) = parse_tsv(&mc);

    // All predefined sequences appear in the output
    assert_eq!(mr.len(), 2);

    // Extra columns appended
    assert!(mh.contains(&"Count".to_owned()));
    assert!(mh.contains(&"Frequency".to_owned()));

    // seq_A: count 3 out of 5 total → frequency 0.6
    let seq_a = mr.iter().find(|r| get(&mh, r, "Name") == "seq_A").unwrap();
    assert_eq!(get(&mh, seq_a, "Count"), "3");
    let freq: f64 = get(&mh, seq_a, "Frequency").parse().unwrap();
    assert!((freq - 0.6).abs() < 1e-5, "frequency should be ~0.6, got {freq}");

    // seq_G: count 1 → frequency 0.2
    let seq_g = mr.iter().find(|r| get(&mh, r, "Name") == "seq_G").unwrap();
    assert_eq!(get(&mh, seq_g, "Count"), "1");

    // Unmatched output
    let uc = fs::read_to_string(&unmatched).unwrap();
    let (uh, ur) = parse_tsv(&uc);

    // CCCC was counted but not in predefined
    assert_eq!(ur.len(), 1);
    assert_eq!(get(&uh, &ur[0], "R1"), "CCCC");
    assert_eq!(get(&uh, &ur[0], "Count"), "1");
}

// Test 7: match — paired-end

#[test]
fn match_paired_basic() {
    let counts_tsv =
        "R1\tR2\tCount\tR1 Name\tR2 Name\nAAAA\tTTTT\t3\tr1\tr1\nCCCC\tGGGG\t1\tr3\tr3\nTTTT\tCCCC\t1\tr5\tr5\n";
    // predefined matches (AAAA,TTTT) and (CCCC,GGGG); (TTTT,CCCC) is unmatched
    let pred_tsv = "Name\tSequence\tSequenceR2\npair_A\tAAAA\tTTTT\npair_C\tCCCC\tGGGG\n";

    let dir = tempfile::tempdir().unwrap();
    let counts = write(dir.path(), "counts.tsv", counts_tsv);
    let predefined = write(dir.path(), "pred.tsv", pred_tsv);
    let matched = dir.path().join("matched.tsv").to_string_lossy().into_owned();
    let unmatched = dir.path().join("unmatched.tsv").to_string_lossy().into_owned();

Command::cargo_bin("mmfqcount").unwrap()
        .args([
            "match",
            "--counts",     &counts,
            "--predefined", &predefined,
            "--seq-col",    "Sequence",
            "--r2-col",     "SequenceR2",
            "--id-col",     "Name",
            "--output",     &matched,
            "--unmatched",  &unmatched,
        ])
        .assert()
        .success();

    let mc = fs::read_to_string(&matched).unwrap();
    let (mh, mr) = parse_tsv(&mc);
    assert_eq!(mr.len(), 2);

    let pair_a = mr.iter().find(|r| get(&mh, r, "Name") == "pair_A").unwrap();
    assert_eq!(get(&mh, pair_a, "Count"), "3");
    let freq: f64 = get(&mh, pair_a, "Frequency").parse().unwrap();
    assert!((freq - 0.6).abs() < 1e-5);

    // (TTTT, CCCC) unmatched
    let uc = fs::read_to_string(&unmatched).unwrap();
    let (uh, ur) = parse_tsv(&uc);
    assert_eq!(ur.len(), 1);
    assert_eq!(get(&uh, &ur[0], "R1"), "TTTT");
    assert_eq!(get(&uh, &ur[0], "R2"), "CCCC");
}

// Test 8: match — predefined sequence with zero count

#[test]
fn match_zero_count_predefined() {
    let counts_tsv = "R1\tCount\tR1 Name\nAAAA\t5\tread1\n";
    // ZZZZ is predefined but never appears in the FASTQ
    let pred_tsv = "Name\tSequence\nseq_A\tAAAA\nseq_Z\tZZZZ\n";

    let dir = tempfile::tempdir().unwrap();
    let counts = write(dir.path(), "counts.tsv", counts_tsv);
    let predefined = write(dir.path(), "pred.tsv", pred_tsv);
    let matched = dir.path().join("matched.tsv").to_string_lossy().into_owned();

Command::cargo_bin("mmfqcount").unwrap()
        .args([
            "match",
            "--counts",     &counts,
            "--predefined", &predefined,
            "--seq-col",    "Sequence",
            "--id-col",     "Name",
            "--output",     &matched,
        ])
        .assert()
        .success();

    let mc = fs::read_to_string(&matched).unwrap();
    let (mh, mr) = parse_tsv(&mc);

    // Both rows present; ZZZZ has count 0 and frequency 0
    let seq_z = mr.iter().find(|r| get(&mh, r, "Name") == "seq_Z").unwrap();
    assert_eq!(get(&mh, seq_z, "Count"), "0");
    let freq: f64 = get(&mh, seq_z, "Frequency").parse().unwrap();
    assert_eq!(freq, 0.0);
}

// Test 9: count — all identical reads

#[test]
fn count_single_all_identical() {
    let fastq = "@r1\nAAAA\n+\nIIII\n@r2\nAAAA\n+\nIIII\n@r3\nAAAA\n+\nIIII\n";

    let dir = tempfile::tempdir().unwrap();
    let r1 = write(dir.path(), "r1.fastq", fastq);
    let out = dir.path().join("out.tsv").to_string_lossy().into_owned();

    Command::cargo_bin("mmfqcount").unwrap()
        .args(["count", "--r1", &r1, "--output", &out])
        .assert()
        .success();

    let content = fs::read_to_string(&out).unwrap();
    let (header, rows) = parse_tsv(&content);

    assert_eq!(rows.len(), 1, "should have exactly one unique sequence");
    assert_eq!(get(&header, &rows[0], "Count"), "3");
    // Example name must be one of the three read names
    let name = get(&header, &rows[0], "R1 Name");
    assert!(["r1", "r2", "r3"].contains(&name));
}

// Test 10: count — gzip input

#[test]
fn count_single_gzip() {
    use std::io::Write as IoWrite;

    let dir = tempfile::tempdir().unwrap();
    let gz_path = dir.path().join("r1.fastq.gz");
    let out = dir.path().join("out.tsv").to_string_lossy().into_owned();

    // Compress SE_FASTQ with flate2
    {
        let f = fs::File::create(&gz_path).unwrap();
        let mut enc = flate2::write::GzEncoder::new(f, flate2::Compression::default());
        enc.write_all(SE_FASTQ.as_bytes()).unwrap();
        enc.finish().unwrap();
    }

Command::cargo_bin("mmfqcount").unwrap()
        .args([
            "count",
            "--r1",    gz_path.to_str().unwrap(),
            "--output", &out,
        ])
        .assert()
        .success();

    let content = fs::read_to_string(&out).unwrap();
    let (header, rows) = parse_tsv(&content);

    // Same result as the plain-text single-end test
    assert_eq!(rows.len(), 3);
    assert_eq!(get(&header, &rows[0], "R1"), "AAAA");
    assert_eq!(get(&header, &rows[0], "Count"), "3");
}

// Test 11: count — single-end, --split-by tag

#[test]
fn count_single_split_by() {
    // Reads carry |sgRNAid=X| tags in read names.
    // AAAA: 2× sgA, 1× sgB  → total 3
    // CCCC: 1× sgB           → total 1
    let fastq = "\
@read1|sgRNAid=sgA\nAAAA\n+\nIIII\n\
@read2|sgRNAid=sgA\nAAAA\n+\nIIII\n\
@read3|sgRNAid=sgB\nAAAA\n+\nIIII\n\
@read4|sgRNAid=sgB\nCCCC\n+\nIIII\n";

    let dir = tempfile::tempdir().unwrap();
    let r1  = write(dir.path(), "r1.fastq", fastq);
    let out = dir.path().join("out.tsv").to_string_lossy().into_owned();

    Command::cargo_bin("mmfqcount").unwrap()
        .args(["count", "--r1", &r1, "--split-by", "sgRNAid", "--output", &out])
        .assert()
        .success();

    let content = fs::read_to_string(&out).unwrap();
    let (header, rows) = parse_tsv(&content);

    // Header should contain per-tag columns (sorted tag values: sgA, sgB)
    assert!(header.contains(&"Count (sgRNAid=sgA)".to_owned()));
    assert!(header.contains(&"Count (sgRNAid=sgB)".to_owned()));
    assert!(header.contains(&"Frequency (sgRNAid=sgA)".to_owned()));
    assert!(header.contains(&"Frequency (sgRNAid=sgB)".to_owned()));

    assert_eq!(rows.len(), 2, "expected 2 unique sequences");

    // AAAA row: 2× sgA, 1× sgB
    let aaaa = rows.iter().find(|r| get(&header, r, "R1") == "AAAA").unwrap();
    assert_eq!(get(&header, aaaa, "Count (sgRNAid=sgA)"), "2");
    assert_eq!(get(&header, aaaa, "Count (sgRNAid=sgB)"), "1");

    // Frequency for sgA: 2 out of 2 sgA reads total → 1.0
    let freq_a: f64 = get(&header, aaaa, "Frequency (sgRNAid=sgA)").parse().unwrap();
    assert!((freq_a - 1.0).abs() < 1e-5, "freq_a={freq_a}");

    // Frequency for sgB: AAAA gets 1 out of 2 sgB reads total → 0.5
    let freq_b: f64 = get(&header, aaaa, "Frequency (sgRNAid=sgB)").parse().unwrap();
    assert!((freq_b - 0.5).abs() < 1e-5, "freq_b={freq_b}");

    // CCCC row: 0× sgA, 1× sgB
    let cccc = rows.iter().find(|r| get(&header, r, "R1") == "CCCC").unwrap();
    assert_eq!(get(&header, cccc, "Count (sgRNAid=sgA)"), "0");
    assert_eq!(get(&header, cccc, "Count (sgRNAid=sgB)"), "1");
}
