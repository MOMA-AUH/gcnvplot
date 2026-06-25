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

Example:

```bash
gcnvplot create-background \
  --read-counts-list background_inputs.txt \
  --output background.tsv

gcnvplot plot \
  --read-counts sample.tsv \
  --background background.tsv \
  --region chr1:100-299 \
  --output plot.svg
```
