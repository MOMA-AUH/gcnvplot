# gcnvplot

[![Anaconda-Server Badge](https://anaconda.org/MOMA-AUH/gcnvplot/badges/version.svg)](https://anaconda.org/MOMA-AUH/gcnvplot) [![Anaconda-Server Badge](https://anaconda.org/MOMA-AUH/gcnvplot/badges/downloads.svg)](https://anaconda.org/MOMA-AUH/gcnvplot)

Tool for plotting GATK germline CNV read-count signals against a background cohort

## Installation

The recommended way to install **gcnvplot** is via [conda](https://docs.conda.io/), using the `MOMA-AUH` channel:

```bash
conda install MOMA-AUH::gcnvplot
```

## Usage

```bash
gcnvplot --help
gcnvplot --version
```

## Inputs and output

`gcnvplot` expects GATK `CollectReadCounts` tables as input. These can be plain TSV files or gzipped TSV files, and must contain the columns `CONTIG`, `START`, `END`, and `COUNT`.

For `create-background`, provide a text file with one sample read-count path per line. This command writes a background cohort TSV with interval-wise normalized summary statistics and a per-interval baseline median.

For `plot`, provide one sample read-count file, one background TSV produced by `create-background`, and a genomic region such as `chr1:100-299`. This command writes an SVG plot showing the sample log2 signal relative to the background cohort.

You can also use `--transcript <TRANSCRIPT_ID>` together with `--gtf <annotations.gtf[.gz]>` to plot by transcript and add an exon track beneath the signal plot. If you want a custom label in the right-side info panel, pass `--sample-name <LABEL>`.

The plot info panel shows `Sample` when provided, then `Gene`, `Transcript`, and `Region`, followed by a separated `Highlight` and `Exons` section when applicable.

Example:

```bash
gcnvplot create-background \
  --read-counts-list background_inputs.txt \
  --output background.tsv

gcnvplot plot \
  --read-counts sample.tsv \
  --background background.tsv \
  --region chr1:100-299 \
  --sample-name SAMPLE_01 \
  --output plot.svg
```

## Synthetic example

A tiny synthetic BRCA1 transcript example is available in [`examples/brca1_synthetic`](examples/brca1_synthetic). It demonstrates a highlighted multi-exon deletion, filled and open sample dots, and an uncovered-exon marker.

![Synthetic BRCA1 transcript plot](examples/brca1_synthetic/brca1_synthetic.svg)

You can render it directly from the repository root:

```bash
gcnvplot plot \
  --read-counts examples/brca1_synthetic/sample_deletion.tsv \
  --background examples/brca1_synthetic/background.tsv \
  --transcript NM_007294.4 \
  --gtf examples/brca1_synthetic/brca1_mane_minimal.gtf \
  --sample-name "Synthetic BRCA1 exon 13-15 deletion" \
  --highlight chr17:43070928-43076614 \
  --output examples/brca1_synthetic/brca1_synthetic.svg
```

This example is designed to show:

- a depressed log2 signal across BRCA1 exons 13-15,
- an intronic interval inside the deletion rendered as an open but still depressed dot,
- an intronic interval outside the deletion rendered as an open near-baseline dot,
- an uncovered exon marked with a triangle in the transcript track.

## Details

`gcnvplot` uses a median-of-ratios normalization. For `create-background`, let `c_ij` be the raw count for interval `i` in background sample `j`.

1. Interval baseline:

   `b_i = median_j(c_ij)` using only positive counts.

2. Background-sample size factor:

   `s_j = median_i(c_ij / b_i)` over intervals with `c_ij > 0` and `b_i > 0`.

3. Normalized background count:

   `n_ij = c_ij / s_j`

The background TSV then stores interval-wise summary statistics across the normalized values `n_ij`, including:

- `BG_NORM_MEAN`
- `BG_NORM_MEDIAN`
- `BG_NORM_SD`
- `BG_NORM_P5`
- `BG_NORM_P95`

For `plot`, let `c_i` be the raw count for the plotted sample at interval `i`.

1. The sample is normalized against the background baselines with the same rule:

   `s = median_i(c_i / b_i)`

   `n_i = c_i / s`

2. The plotted signal is the stabilized log2 ratio against the background median:

   `signal_i = log2((n_i + 0.01) / (m_i + 0.01))`

   where `m_i = BG_NORM_MEDIAN`.

3. The background ribbon is drawn by transforming the stored background percentiles in the same way:

   `lower_i = log2((p5_i + 0.01) / (m_i + 0.01))`

   `upper_i = log2((p95_i + 0.01) / (m_i + 0.01))`

This means the plotted curve shows relative dosage after library-size normalization, while the ribbon shows where the central background cohort typically lies for each interval.
