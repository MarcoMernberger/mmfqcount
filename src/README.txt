Benutzung


Paired-end

fastq_counter count \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  --output counts.tsv


Mit Trimming:

fastq_counter count \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  --trim-start ATGCGT \
  --trim-stop  TTAGCA \
  --trim-length 20 \
  --output counts.tsv


Single-end

fastq_counter count \
  -1 sample_R1.fastq.gz \
  --output counts.tsv


Predefined sequences matchen

fastq_counter match \
  --counts counts.tsv \
  --predefined sequences.tsv \
  --seq-col Sequence \
  --id-col Name \
  --output matched.tsv \
  --unmatched unmatched.tsv


Für paired matching mit R2-Spalte in der predefined-Datei:

fastq_counter match \
  --counts counts.tsv \
  --predefined sequences.tsv \
  --seq-col Sequence \
  --r2-col Sequence_R2 \
  --id-col Name \
  --output matched.tsv \
  --unmatched unmatched.tsv
