"""Self-contained HTML report generation."""

from __future__ import annotations

import csv
import html
import json
import os
import tempfile
from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from .models import TranscriptAnnotation
from .plotting import prepare_plot_data, render_prepared_plot_svg
from .segments import CnvSegment, parse_segments_vcf, segments_for_transcript
from .transcripts import TranscriptIndex, normalize_gene_identifier


@dataclass(frozen=True)
class GeneRequest:
    """One requested gene and its optional explicit transcript."""

    gene_name: str
    stable_id: str
    transcript_id: str | None = None


@dataclass(frozen=True)
class ReportBuildResult:
    """Summary of a successfully built report."""

    gene_count: int
    segment_count: int
    vcf_sample: str | None


def read_gene_requests(path: Path) -> list[GeneRequest]:
    """Read and validate a tab-separated report gene list."""
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"gene_name", "stable_id"}
        fields = set(reader.fieldnames or [])
        missing = required - fields
        if missing:
            raise ValueError(
                f"{path}: missing required column(s): {', '.join(sorted(missing))}"
            )

        requests: list[GeneRequest] = []
        seen_names: set[str] = set()
        seen_ids: set[str] = set()
        for line_number, row in enumerate(reader, start=2):
            gene_name = (row.get("gene_name") or "").strip()
            if not gene_name:
                raise ValueError(f"{path}:{line_number}: gene_name is empty")
            stable_id = normalize_gene_identifier(row.get("stable_id") or "")
            if not stable_id:
                raise ValueError(f"{path}:{line_number}: stable_id is empty")
            name_key = gene_name.casefold()
            id_key = stable_id.casefold()
            if name_key in seen_names:
                raise ValueError(
                    f"{path}:{line_number}: duplicate gene name {gene_name!r}"
                )
            if id_key in seen_ids:
                raise ValueError(
                    f"{path}:{line_number}: duplicate stable identifier {stable_id!r}"
                )
            transcript = (row.get("transcript") or "").strip() or None
            requests.append(GeneRequest(gene_name, stable_id, transcript))
            seen_names.add(name_key)
            seen_ids.add(id_key)

    if not requests:
        raise ValueError(f"{path}: no genes found")
    return requests


def _sample_label(path: Path) -> str:
    name = path.name
    for suffix in (".tsv.gz", ".txt.gz", ".tsv", ".txt", ".gz"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def _tool_version() -> str:
    try:
        return version("gcnvplot")
    except PackageNotFoundError:
        return "development"


def _report_html(
    plots: list[tuple[GeneRequest, str, str, list[CnvSegment]]],
    *,
    sample_name: str,
    annotation_build: str | None,
) -> str:
    templates: list[str] = []
    lookup: dict[str, object] = {}
    segment_count = 0
    for number, (request, transcript_id, svg, segments) in enumerate(plots):
        template_id = f"plot-{number}"
        segment_count += len(segments)
        lookup[request.gene_name.upper()] = {
            "template": template_id,
            "gene": request.gene_name,
            "stableId": request.stable_id,
            "transcript": transcript_id,
            "segments": [
                {
                    "label": segment.label,
                    "type": segment.svtype,
                    "filtered": segment.filter_status not in {".", "PASS"},
                }
                for segment in segments
            ],
        }
        templates.append(f'<template id="{template_id}">{svg}</template>')

    lookup_json = json.dumps(lookup, ensure_ascii=True)
    for character, escaped in (("<", "\\u003c"), (">", "\\u003e"), ("&", "\\u0026")):
        lookup_json = lookup_json.replace(character, escaped)
    provenance = [f"gcnvplot {_tool_version()}"]
    if annotation_build is not None:
        provenance.append(annotation_build)
    provenance_text = " · ".join(provenance)

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta http-equiv="Content-Security-Policy"
    content="default-src 'none'; script-src 'unsafe-inline'; style-src 'unsafe-inline'">
  <title>{html.escape(sample_name)} — gcnvplot report</title>
  <style>
    :root {{ color-scheme: light; font-family: Inter, system-ui, sans-serif; }}
    body {{ margin: 0; background: #f3f6f8; color: #17212b; }}
    main {{ max-width: 1480px; margin: 0 auto; padding: 2.5rem 1.5rem; }}
    h1 {{ margin: 0 0 .4rem; font-size: 1.8rem; }}
    #details, .provenance {{ color: #536170; }}
    .provenance {{ font-size: .82rem; margin-top: .5rem; }}
    form {{ display: flex; gap: .6rem; max-width: 38rem; margin: 1.8rem 0 .7rem; position: relative; }}
    label {{ position: absolute; width: 1px; height: 1px; overflow: hidden; }}
    input {{ flex: 1; min-width: 0; padding: .75rem .9rem; border: 1px solid #aeb8c2;
      border-radius: .4rem; font: inherit; text-transform: uppercase; }}
    input:focus {{ outline: 3px solid #b9dcf4; border-color: #246b9c; }}
    button {{ padding: .75rem 1.1rem; border: 0; border-radius: .4rem;
      background: #126a8a; color: white; font: inherit; font-weight: 700; cursor: pointer; }}
    #gene-suggestions {{ display: none; position: absolute; z-index: 2; top: calc(100% + .2rem);
      left: 0; right: 6.8rem; max-height: 14rem; overflow-y: auto; list-style: none;
      margin: 0; padding: .25rem; border: 1px solid #aeb8c2; border-radius: .4rem;
      background: white; box-shadow: 0 4px 12px rgb(0 0 0 / 16%); }}
    #gene-suggestions.visible {{ display: block; }}
    #gene-suggestions li {{ padding: .55rem .65rem; border-radius: .25rem; cursor: pointer; }}
    #gene-suggestions li:hover, #gene-suggestions li[aria-selected="true"] {{
      background: #dceef8; }}
    #message {{ min-height: 1.5rem; margin: 0 0 1rem; color: #a11b1b; font-weight: 650; }}
    #details {{ display: none; margin: 0 0 .8rem; }}
    #variants {{ display: none; gap: .45rem; flex-wrap: wrap; list-style: none;
      padding: 0; margin: 0 0 1rem; }}
    #variants li {{ padding: .3rem .55rem; border-radius: 999px; font-size: .85rem;
      border: 1px solid; }}
    #variants .DEL {{ color: #991b1b; background: #fee2e2; border-color: #fca5a5; }}
    #variants .DUP {{ color: #1e40af; background: #dbeafe; border-color: #93c5fd; }}
    #variants .filtered {{ border-style: dashed; }}
    #plot {{ display: none; background: white; padding: 1rem; border-radius: .55rem;
      box-shadow: 0 1px 5px rgb(0 0 0 / 12%); overflow-x: auto; }}
    #plot svg {{ display: block; width: 100%; height: auto; min-width: 900px; }}
    @media print {{
      body {{ background: white; }}
      form, #message, .provenance {{ display: none; }}
      #plot {{ box-shadow: none; padding: 0; overflow: visible; }}
      #plot svg {{ min-width: 0; }}
    }}
  </style>
</head>
<body>
<main>
  <header>
    <h1>{html.escape(sample_name)}</h1>
    <p class="provenance">{html.escape(provenance_text)}</p>
  </header>
  <noscript>This report requires JavaScript to select an embedded plot.</noscript>
  <form id="search-form" novalidate>
    <label for="gene-query">Gene name</label>
    <input id="gene-query" name="report-gene-query" type="text" autocomplete="off"
      autocapitalize="characters" spellcheck="false" placeholder="Gene name" required autofocus
      role="combobox" aria-autocomplete="list" aria-controls="gene-suggestions" aria-expanded="false">
    <button type="submit">Show plot</button>
    <ul id="gene-suggestions" role="listbox" aria-label="Matching report genes"></ul>
  </form>
  <p id="message" role="alert" aria-live="polite"></p>
  <p id="details"></p>
  <ul id="variants" aria-label="Overlapping CNV segments"></ul>
  <section id="plot" aria-live="polite"></section>
  {"".join(templates)}
</main>
<script>
  const genes = {lookup_json};
  const form = document.getElementById("search-form");
  const query = document.getElementById("gene-query");
  const message = document.getElementById("message");
  const details = document.getElementById("details");
  const variants = document.getElementById("variants");
  const plot = document.getElementById("plot");
  const suggestions = document.getElementById("gene-suggestions");
  const geneNames = Object.values(genes).map((entry) => entry.gene)
    .sort((left, right) => left.localeCompare(right));
  let activeSuggestion = -1;

  function hideSuggestions() {{
    suggestions.replaceChildren();
    suggestions.classList.remove("visible");
    query.setAttribute("aria-expanded", "false");
    query.removeAttribute("aria-activedescendant");
    activeSuggestion = -1;
  }}

  function matchingGenes(value) {{
    const typed = value.trim().toUpperCase();
    if (!typed) return [];
    const startsWith = geneNames.filter((gene) => gene.toUpperCase().startsWith(typed));
    const contains = geneNames.filter((gene) =>
      !gene.toUpperCase().startsWith(typed) && gene.toUpperCase().includes(typed));
    return startsWith.concat(contains).slice(0, 8);
  }}

  function chooseSuggestion(gene) {{
    query.value = gene;
    hideSuggestions();
    showGene(gene);
    query.focus();
  }}

  function renderSuggestions() {{
    const matches = matchingGenes(query.value);
    suggestions.replaceChildren();
    activeSuggestion = -1;
    if (!matches.length) {{
      suggestions.classList.remove("visible");
      query.setAttribute("aria-expanded", "false");
      return;
    }}
    matches.forEach((gene, index) => {{
      const item = document.createElement("li");
      item.id = `gene-suggestion-${{index}}`;
      item.setAttribute("role", "option");
      item.setAttribute("aria-selected", "false");
      item.textContent = gene;
      item.addEventListener("mousedown", (event) => {{
        event.preventDefault();
        chooseSuggestion(gene);
      }});
      suggestions.appendChild(item);
    }});
    suggestions.classList.add("visible");
    query.setAttribute("aria-expanded", "true");
  }}

  function setActiveSuggestion(index) {{
    const items = suggestions.children;
    if (!items.length) return;
    activeSuggestion = (index + items.length) % items.length;
    Array.from(items).forEach((item, itemIndex) => {{
      item.setAttribute("aria-selected", String(itemIndex === activeSuggestion));
    }});
    const active = items[activeSuggestion];
    query.setAttribute("aria-activedescendant", active.id);
    active.scrollIntoView({{ block: "nearest" }});
  }}

  function showGene(value) {{
    const typed = value.trim();
    const entry = genes[typed.toUpperCase()];
    plot.replaceChildren();
    plot.style.display = "none";
    details.style.display = "none";
    variants.replaceChildren();
    variants.style.display = "none";
    if (!entry) {{
      message.textContent = typed
        ? `Gene "${{typed}}" is not included in this report.`
        : "Type a gene name.";
      return;
    }}
    query.value = entry.gene;
    message.textContent = "";
    details.textContent =
      `${{entry.gene}} · ${{entry.stableId}} · ${{entry.transcript}}` +
      (entry.segments.length ? ` · ${{entry.segments.length}} CNV segment(s)` : "");
    details.style.display = "block";
    for (const segment of entry.segments) {{
      const item = document.createElement("li");
      item.className = segment.type + (segment.filtered ? " filtered" : "");
      item.textContent = segment.label;
      variants.appendChild(item);
    }}
    if (entry.segments.length) variants.style.display = "flex";
    const template = document.getElementById(entry.template);
    plot.appendChild(template.content.cloneNode(true));
    plot.style.display = "block";
    history.replaceState(null, "", "#" + encodeURIComponent(entry.gene));
  }}

  form.addEventListener("submit", (event) => {{
    event.preventDefault();
    hideSuggestions();
    showGene(query.value);
  }});
  query.addEventListener("input", renderSuggestions);
  query.addEventListener("focus", renderSuggestions);
  query.addEventListener("blur", () => window.setTimeout(hideSuggestions, 100));
  query.addEventListener("keydown", (event) => {{
    if (event.key === "ArrowDown") {{
      if (!suggestions.children.length) renderSuggestions();
      if (suggestions.children.length) {{
        event.preventDefault();
        setActiveSuggestion(activeSuggestion + 1);
      }}
    }} else if (event.key === "ArrowUp") {{
      if (suggestions.children.length) {{
        event.preventDefault();
        setActiveSuggestion(activeSuggestion - 1);
      }}
    }} else if (event.key === "Enter" && activeSuggestion >= 0) {{
      event.preventDefault();
      chooseSuggestion(suggestions.children[activeSuggestion].textContent);
    }} else if (event.key === "Escape") {{
      hideSuggestions();
    }}
  }});
  if (location.hash.length > 1) showGene(decodeURIComponent(location.hash.slice(1)));
</script>
</body>
</html>
"""


def _atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
        text=True,
    )
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            handle.write(text)
        os.replace(temporary_path, path)
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise


def build_report(
    read_counts_path: Path,
    background_path: Path,
    genes_path: Path,
    transcript_db_path: Path,
    output_path: Path,
    *,
    sample_name: str | None = None,
    strict_intervals: bool = True,
    segments_vcf_path: Path | None = None,
    include_filtered_segments: bool = False,
) -> ReportBuildResult:
    """Build one searchable, self-contained HTML gene report."""
    requests = read_gene_requests(genes_path)
    prepared = prepare_plot_data(read_counts_path, background_path)
    shared_contigs = set(prepared.counts_by_contig) & set(prepared.background_by_contig)
    if not shared_contigs:
        raise ValueError("read counts and background do not share any contig names")
    selected_vcf_sample: str | None = None
    all_segments: list[CnvSegment] = []
    if segments_vcf_path is not None:
        selected_vcf_sample, all_segments = parse_segments_vcf(
            segments_vcf_path,
            sample_name=sample_name,
            include_filtered=include_filtered_segments,
        )

    label = sample_name or selected_vcf_sample or _sample_label(read_counts_path)
    resolved: list[tuple[GeneRequest, TranscriptAnnotation]] = []
    errors: list[str] = []
    annotation_build: str | None = None
    with TranscriptIndex(transcript_db_path) as transcript_index:
        index_metadata = transcript_index.metadata()
        annotation_build = index_metadata.get("genome-build")
        for request in requests:
            try:
                annotation = transcript_index.resolve_gene(
                    request.gene_name,
                    request.stable_id,
                    request.transcript_id,
                    available_contigs=shared_contigs,
                )
                resolved.append((request, annotation))
            except ValueError as error:
                errors.append(str(error))
        if errors:
            raise ValueError("Report gene validation failed:\n- " + "\n- ".join(errors))

        plots: list[tuple[GeneRequest, str, str, list[CnvSegment]]] = []
        for request, annotation in resolved:
            transcript_segments = segments_for_transcript(all_segments, annotation)
            overlays = [
                segment.as_overlay(contig=annotation.contig)
                for segment in transcript_segments
            ]
            try:
                svg = render_prepared_plot_svg(
                    prepared,
                    (annotation.contig, annotation.start, annotation.end),
                    transcript=annotation,
                    sample_name=label,
                    overlays=overlays,
                    strict_intervals=strict_intervals,
                )
            except ValueError as error:
                errors.append(f"{request.gene_name}: {error}")
                continue
            plots.append(
                (
                    request,
                    annotation.transcript_id,
                    svg,
                    transcript_segments,
                )
            )
    if errors:
        raise ValueError("Report plot validation failed:\n- " + "\n- ".join(errors))

    page = _report_html(
        plots,
        sample_name=label,
        annotation_build=annotation_build,
    )
    _atomic_write_text(output_path, page)
    return ReportBuildResult(
        gene_count=len(plots),
        segment_count=sum(len(item[3]) for item in plots),
        vcf_sample=selected_vcf_sample,
    )
