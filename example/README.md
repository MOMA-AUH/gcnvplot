# BRCA1 synthetic plotting example

This example is a tiny, synthetic dataset for exercising the transcript plotting features in `gcnvplot`.

It contains:

- `brca1_mane_minimal.gtf`: the BRCA1 `NM_007294.4` transcript line plus exon lines from an hg38 NCBI RefSeq GTF.
- `background.tsv`: synthetic interval-wise background statistics with 20 mock background samples and interval-to-interval depth variation.
- `sample_deletion.tsv`: a synthetic sample with realistic count jitter and a heterozygous deletion across BRCA1 exons 13-15.

The intervals are roughly 1 kb and have different expected depths, mimicking variable probe efficiency. Two intervals are intentionally intronic so the plot shows open sample dots; one of those intronic intervals falls inside the deletion and is depressed with the deleted exons. Exon 12 has no overlapping interval so the transcript track shows the uncovered-exon triangle.

Run from the repository root:

```bash
gcnvplot index-gtf \
  --gtf example/brca1_mane_minimal.gtf \
  --output example/brca1_mane_minimal.sqlite

gcnvplot plot \
  --read-counts example/sample_deletion.tsv \
  --background example/background.tsv \
  --transcript NM_007294.4 \
  --transcript-db example/brca1_mane_minimal.sqlite \
  --sample-name "Synthetic BRCA1 exon 13-15 deletion" \
  --highlight chr17:43070928-43076614 \
  --output example/brca1_synthetic.svg
```

The output should show:

- a depressed log2 signal across the highlighted exon 13-15 deletion, including the intronic interval inside the deletion,
- filled dots for intervals overlapping exons,
- open dots for intervals outside exons,
- a red triangle under exon 12, which has no overlapping synthetic interval.
