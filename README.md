# gatk-germline-cnv-plotter

Basic scaffold for a Python CLI application for plotting GATK germline CNV outputs.

## Development

Install the project in editable mode with test dependencies:

```bash
python -m pip install -e .[dev]
```

Run the test suite:

```bash
pytest
```

Try the CLI:

```bash
gatk-germline-cnv-plotter --version
gatk-germline-cnv-plotter plot --input sample.tsv --output plot.png
```
