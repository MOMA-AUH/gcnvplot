# gcnvplot

Basic scaffold for a Python CLI application for plotting GATK germline CNV outputs.

## Development

Install the project in editable mode with test dependencies:

```bash
python -m pip install -e .[dev]
```

Or create the reference conda development environment:

```bash
conda env create -f environment.yml
conda activate gcnvplot-dev
```

Run the test suite:

```bash
pytest
```

Try the CLI:

```bash
gcnvplot --version
gcnvplot plot --input sample.tsv --output plot.png
```
