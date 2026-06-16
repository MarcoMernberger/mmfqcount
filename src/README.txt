Usage


Paired-end

mmfqcount count \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  --output counts.tsv


With Trimming:

mmfqcount count \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  --trim-start ATGCGT \
  --trim-stop  TTAGCA \
  --trim-length 20 \
  --output counts.tsv


Single-end

mmfqcount count \
  -1 sample_R1.fastq.gz \
  --output counts.tsv


Predefined sequences matching

mmfqcount match \
  --counts counts.tsv \
  --predefined sequences.tsv \
  --seq-col Sequence \
  --id-col Name \
  --output matched.tsv \
  --unmatched unmatched.tsv


 paired matching with R2-Spalte in  predefined-file:

mmfqcount match \
  --counts counts.tsv \
  --predefined sequences.tsv \
  --seq-col Sequence \
  --r2-col Sequence_R2 \
  --id-col Name \
  --output matched.tsv \
  --unmatched unmatched.tsv
