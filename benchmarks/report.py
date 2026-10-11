"""Render benchmark results as Markdown and HTML.

The report is the deliverable a reader actually looks
at, so it leads with the comparison table and then
breaks each metric out. Numbers are formatted to a
fixed precision and ``None`` (a metric that could not
be computed, e.g. specificity without a negative set)
is rendered as ``--`` rather than as a misleading 0.

Markdown is the default because it drops straight into
a PR or an issue; HTML is offered for the standalone
dashboard. Both are generated from the same
:class:`BenchmarkResult` list, so they cannot disagree.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from .amr_benchmark import BenchmarkResult

__all__ = ["render_markdown", "render_html", "write_report"]


def _fmt(value: float | None, precision: int = 3) -> str:
    if value is None:
        return "--"
    if isinstance(value, float) and (value != value):  # NaN
        return "--"
    return f"{value:.{precision}f}"


def _pct(value: float | None) -> str:
    if value is None:
        return "--"
    return f"{value * 100:.1f}%"


def render_markdown(results: Sequence[BenchmarkResult]) -> str:
    """A comparison table plus per-metric breakdown."""
    lines: list[str] = []
    lines.append("# AMR benchmark results\n")

    lines.append("## Summary\n")
    lines.append(
        "| tool | version | sensitivity | precision | F1 | "
        "epistasis F1 | dN/dS r | lineage | runtime | "
        "genomes/hr |"
    )
    lines.append(
        "|------|---------|-------------|-----------|-----|"
        "--------------|---------|---------|---------|"
        "-----------|"
    )
    for result in results:
        lines.append(
            f"| {result.tool} "
            f"| {result.version or '--'} "
            f"| {_pct(result.mutations.sensitivity)} "
            f"| {_pct(result.mutations.precision)} "
            f"| {_fmt(result.mutations.f1)} "
            f"| {_fmt(result.epistasis.f1)} "
            f"| {_fmt(result.dn_ds.pearson)} "
            f"| {_pct(result.lineage.concordance)} "
            f"| {result.performance.runtime_seconds:.1f}s "
            f"| {result.performance.genomes_per_hour:.0f} |"
        )
    lines.append("")

    lines.append("## Mutation calling\n")
    lines.append("| tool | TP | FP | FN | TN | novel recall |")
    lines.append("|------|----|----|----|----|--------------|")
    for result in results:
        m = result.mutations
        lines.append(
            f"| {result.tool} | {m.tp} | {m.fp} | {m.fn} | "
            f"{m.tn or '--'} | {_pct(m.novel_recall)} |"
        )
    lines.append("")

    lines.append("## Epistasis\n")
    lines.append("| tool | TP | FP | FN | precision | recall | F1 |")
    lines.append("|------|----|----|----|-----------|--------|-----|")
    for result in results:
        e = result.epistasis
        lines.append(
            f"| {result.tool} | {e.tp} | {e.fp} | {e.fn} | "
            f"{_pct(e.precision)} | {_pct(e.recall)} | "
            f"{_fmt(e.f1)} |"
        )
    lines.append("")

    lines.append("## dN/dS\n")
    lines.append("| tool | Pearson | Spearman | MAE | genes |")
    lines.append("|------|---------|----------|-----|-------|")
    for result in results:
        d = result.dn_ds
        lines.append(
            f"| {result.tool} | {_fmt(d.pearson)} | "
            f"{_fmt(d.spearman)} | {_fmt(d.mae)} | {d.n_genes} |"
        )
    lines.append("")

    lines.append("## Performance\n")
    lines.append("| tool | runtime | peak RSS | genomes | genomes/hr |")
    lines.append("|------|---------|----------|---------|------------|")
    for result in results:
        p = result.performance
        memory = f"{p.peak_memory_mb:.0f} MiB" if p.peak_memory_mb else "--"
        lines.append(
            f"| {result.tool} | {p.runtime_seconds:.1f}s | {memory} | "
            f"{p.n_genomes} | {p.genomes_per_hour:.0f} |"
        )
    lines.append("")
    return "\n".join(lines)


def render_html(results: Sequence[BenchmarkResult]) -> str:
    """A standalone HTML dashboard of the same numbers."""
    rows = []
    for result in results:
        rows.append(_html_row(result))
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>AMR benchmark results</title>
<style>
  body {{ font-family: system-ui, sans-serif; margin: 2em; }}
  table {{ border-collapse: collapse; margin: 1em 0; }}
  th, td {{ border: 1px solid #ccc; padding: 0.4em 0.8em;
            text-align: right; }}
  th {{ background: #f4f4f4; }}
  td:first-child, th:first-child {{ text-align: left; }}
</style>
</head>
<body>
<h1>AMR benchmark results</h1>
<table>
<thead><tr>
  <th>tool</th><th>sensitivity</th><th>precision</th>
  <th>F1</th><th>epistasis F1</th><th>dN/dS r</th>
  <th>lineage</th><th>runtime</th><th>genomes/hr</th>
</tr></thead>
<tbody>
{chr(10).join(rows)}
</tbody>
</table>
</body>
</html>
"""


def _html_row(result: BenchmarkResult) -> str:
    m = result.mutations
    p = result.performance
    return (
        "<tr>"
        f"<td>{result.tool}</td>"
        f"<td>{_pct(m.sensitivity)}</td>"
        f"<td>{_pct(m.precision)}</td>"
        f"<td>{_fmt(m.f1)}</td>"
        f"<td>{_fmt(result.epistasis.f1)}</td>"
        f"<td>{_fmt(result.dn_ds.pearson)}</td>"
        f"<td>{_pct(result.lineage.concordance)}</td>"
        f"<td>{p.runtime_seconds:.1f}s</td>"
        f"<td>{p.genomes_per_hour:.0f}</td>"
        "</tr>"
    )


def write_report(
    results: Sequence[BenchmarkResult],
    output_dir: Any,
) -> dict[str, str]:
    """Write both formats and return their paths."""
    from pathlib import Path

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    markdown = render_markdown(results)
    html = render_html(results)
    md_path = output_dir / "benchmark_report.md"
    html_path = output_dir / "benchmark_report.html"
    md_path.write_text(markdown, encoding="utf-8")
    html_path.write_text(html, encoding="utf-8")
    return {"markdown": str(md_path), "html": str(html_path)}
