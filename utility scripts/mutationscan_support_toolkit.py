import argparse
import csv
import json
import re
import shutil
import time
from collections import Counter
from pathlib import Path

import pandas as pd
import requests

API_GENOME = "https://www.bv-brc.org/api/genome/"
API_GENOME_SEQ = "https://www.bv-brc.org/api/genome_sequence/"


def _normalize_ids(series):
    seen = set()
    ordered = []
    for raw in series:
        gid = str(raw).strip()
        if not gid or gid.lower() == "nan":
            continue
        if gid in seen:
            continue
        seen.add(gid)
        ordered.append(gid)
    return ordered


def _request_json(url, timeout, retries):
    last_err = "unknown"
    for attempt in range(1, retries + 1):
        try:
            r = requests.get(url, timeout=timeout)
            if r.status_code == 200:
                return r.json(), None
            last_err = f"HTTP {r.status_code}"
        except requests.RequestException as exc:
            last_err = f"{type(exc).__name__}: {exc}"
        time.sleep(0.25 * attempt)
    return None, last_err


def cmd_download_rest(args):
    output_dir = Path(args.output_dir)
    csv_file = Path(args.csv_file)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(csv_file, dtype=str, keep_default_na=False)
    if args.genome_id_column not in df.columns:
        raise KeyError(f"Column '{args.genome_id_column}' missing in {csv_file}")

    all_ids = _normalize_ids(df[args.genome_id_column])

    downloaded_ids = set()
    for f in output_dir.glob("*.fna"):
        try:
            if f.stat().st_size > args.min_bytes:
                downloaded_ids.add(f.stem)
        except OSError:
            continue

    missing = [gid for gid in all_ids if gid not in downloaded_ids]
    print(f"Total IDs in CSV: {len(all_ids)}")
    print(f"Missing IDs before REST download: {len(missing)}")

    if not missing:
        print("All genomes already present. Nothing to download.")
        return

    headers = {"Accept": "application/dna+fasta", "User-Agent": "MutationScan-Toolkit/1.0"}

    def api_url(gid):
        return f"{API_GENOME_SEQ}?eq(genome_id,{gid})&limit({args.api_limit})"

    success = 0
    failed = []

    for i, gid in enumerate(missing, start=1):
        out_path = output_dir / f"{gid}.fna"
        if out_path.exists() and out_path.stat().st_size > args.min_bytes:
            continue

        last_err = "unknown"
        ok = False
        for attempt in range(1, args.retries + 1):
            try:
                r = requests.get(api_url(gid), headers=headers, timeout=args.timeout)
                if r.status_code != 200:
                    last_err = f"HTTP {r.status_code}"
                    time.sleep(0.2 * attempt)
                    continue
                payload = r.text or ""
                if len(payload.encode("utf-8", errors="ignore")) <= args.min_bytes:
                    last_err = "payload_too_small"
                    time.sleep(0.2 * attempt)
                    continue
                tmp = output_dir / f"{gid}.fna.part"
                tmp.write_text(payload, encoding="utf-8")
                tmp.replace(out_path)
                success += 1
                ok = True
                break
            except requests.RequestException as exc:
                last_err = f"{type(exc).__name__}: {exc}"
                time.sleep(0.2 * attempt)
        if not ok:
            failed.append((gid, last_err))

        if i % 25 == 0 or i == len(missing):
            print(f"Progress {i}/{len(missing)} | success={success} failed={len(failed)}")

    failed_log = output_dir / args.failed_log
    missing_log = output_dir / args.missing_log
    failed_log.write_text("\n".join(g for g, _ in failed) + ("\n" if failed else ""), encoding="utf-8")

    present = {f.stem for f in output_dir.glob("*.fna") if f.stat().st_size > args.min_bytes}
    remaining = [gid for gid in all_ids if gid not in present]
    missing_log.write_text("\n".join(remaining) + ("\n" if remaining else ""), encoding="utf-8")

    print("=" * 60)
    print(f"Newly downloaded: {success}")
    print(f"Failed this run: {len(failed)}")
    print(f"Remaining missing after run: {len(remaining)}")
    print(f"Failed IDs log: {failed_log}")
    print(f"Remaining missing log: {missing_log}")
    print("=" * 60)


def cmd_fetch_metadata(args):
    input_csv = Path(args.input_csv)
    df = pd.read_csv(input_csv, dtype=str, keep_default_na=False)
    if args.genome_id_column not in df.columns:
        raise KeyError(f"Column '{args.genome_id_column}' missing in {input_csv}")

    out_meta = Path(args.output_metadata_csv)
    out_enriched = Path(args.output_enriched_csv)
    failed_log = Path(args.failed_log)
    out_meta.parent.mkdir(parents=True, exist_ok=True)
    out_enriched.parent.mkdir(parents=True, exist_ok=True)
    failed_log.parent.mkdir(parents=True, exist_ok=True)

    ids = _normalize_ids(df[args.genome_id_column])

    def batch_query(batch):
        quoted = ",".join([f'"{gid}"' for gid in batch])
        limit = max(200, len(batch) * 3)
        url = f"{API_GENOME}?in(genome_id,({quoted}))&limit({limit})"
        payload, err = _request_json(url, timeout=args.timeout, retries=args.retries)
        if payload is None:
            return {}, err
        out = {}
        for item in payload:
            gid = str(item.get("genome_id", "")).strip()
            if gid:
                out[gid] = item
        return out, None

    def single_query(gid):
        url = f"{API_GENOME}?eq(genome_id,{gid})&limit(1)"
        payload, err = _request_json(url, timeout=args.timeout, retries=args.retries)
        if payload is None:
            return None, err
        if not payload:
            return None, "not_found"
        return payload[0], None

    meta = {}
    failures = []
    batches = [ids[i : i + args.batch_size] for i in range(0, len(ids), args.batch_size)]

    for i, batch in enumerate(batches, start=1):
        result, err = batch_query(batch)
        if err is not None:
            result = {}
        for gid in batch:
            if gid in result:
                meta[gid] = result[gid]
            else:
                item, e = single_query(gid)
                if item is not None:
                    meta[gid] = item
                else:
                    failures.append((gid, e or "unknown"))
        if i % 10 == 0 or i == len(batches):
            print(f"Progress {i}/{len(batches)} batches | metadata_found={len(meta)} failures={len(failures)}")
        time.sleep(args.sleep_seconds)

    keys = {"genome_id", "metadata_found"}
    for item in meta.values():
        keys.update(item.keys())
    cols = ["genome_id", "metadata_found"] + sorted(k for k in keys if k not in {"genome_id", "metadata_found"})

    rows = []
    for gid in ids:
        if gid in meta:
            row = {"genome_id": gid, "metadata_found": 1}
            for k, v in meta[gid].items():
                if isinstance(v, (dict, list)):
                    row[k] = json.dumps(v, ensure_ascii=True)
                else:
                    row[k] = v
        else:
            row = {"genome_id": gid, "metadata_found": 0}
        rows.append(row)

    with out_meta.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    meta_df = pd.DataFrame(rows).rename(columns={"genome_id": args.genome_id_column})
    enriched = df.merge(meta_df, on=args.genome_id_column, how="left")
    enriched.to_csv(out_enriched, index=False)

    failed_log.write_text("\n".join(f"{gid}\t{err}" for gid, err in failures) + ("\n" if failures else ""), encoding="utf-8")

    print("=" * 60)
    print(f"Metadata found: {len(meta)} / {len(ids)}")
    print(f"Failures: {len(failures)}")
    print(f"Metadata CSV: {out_meta}")
    print(f"Enriched CSV: {out_enriched}")
    print(f"Failed IDs log: {failed_log}")
    print("=" * 60)


def _clean_location(value):
    if pd.isna(value):
        return "Unknown"
    text = str(value).strip()
    if not text:
        return "Unknown"
    if ":" in text:
        text = text.split(":", 1)[1].strip()
    text = re.sub(r"\s*\([^)]*\)", "", text)
    text = re.sub(r"\s+", " ", text).strip(" ,;-")
    return text if text else "Unknown"


def cmd_geospatial_matrix(args):
    keep_cols = [
        "Genome ID",
        "Antibiotic",
        "Resistant Phenotype",
        "isolation_country",
        "geographic_location",
        "latitude",
        "longitude",
        "collection_year",
        "host_group",
        "isolation_source",
        "disease",
    ]

    frames = []
    for p in args.metadata_csvs:
        path = Path(p)
        if not path.exists():
            raise FileNotFoundError(f"Metadata CSV not found: {path}")
        df = pd.read_csv(path, dtype=str, keep_default_na=False)
        for c in keep_cols:
            if c not in df.columns:
                df[c] = ""
        frames.append(df[keep_cols].copy())

    meta = pd.concat(frames, ignore_index=True)
    meta["Genome ID"] = meta["Genome ID"].astype(str).str.strip()
    meta = meta[meta["Genome ID"] != ""].drop_duplicates(subset=["Genome ID"], keep="first").reset_index(drop=True)
    meta["geographic_location_clean"] = meta["geographic_location"].apply(_clean_location)

    rg = {g.strip().lower() for g in args.regulatory_genes if g.strip()}
    g_frames = []
    for p in args.genomics_reports:
        path = Path(p)
        if not path.exists():
            raise FileNotFoundError(f"Genomics report not found: {path}")
        gdf = pd.read_csv(path, dtype=str, keep_default_na=False)
        for c in ["Accession", "Gene", "Mutation"]:
            if c not in gdf.columns:
                raise KeyError(f"Missing {c} in {path}")
        gdf["Accession"] = gdf["Accession"].astype(str).str.strip()
        gdf["Gene"] = gdf["Gene"].astype(str).str.strip()
        gdf["Mutation"] = gdf["Mutation"].astype(str).str.strip()
        gdf = gdf[(gdf["Accession"] != "") & (gdf["Gene"] != "") & (gdf["Mutation"] != "")]
        gdf = gdf[gdf["Gene"].str.lower().isin(rg)]
        g_frames.append(gdf)

    genomics = pd.concat(g_frames, ignore_index=True).drop_duplicates(subset=["Accession", "Gene", "Mutation"]).reset_index(drop=True)

    merged = meta.merge(genomics, left_on="Genome ID", right_on="Accession", how="inner")

    freq = (
        merged.groupby(["geographic_location_clean", "Mutation"], dropna=False)["Genome ID"]
        .nunique()
        .reset_index(name="mutation_frequency")
    )
    matrix = freq.pivot(index="geographic_location_clean", columns="Mutation", values="mutation_frequency").fillna(0).astype(int)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    meta.to_csv(out_dir / "Master_Geospatial_AMR.csv", index=False)
    merged.to_csv(out_dir / "Merged_Regulatory_Geospatial_AMR.csv", index=False)
    freq.to_csv(out_dir / "Geospatial_Mutation_Long.csv", index=False)
    matrix.to_csv(out_dir / "Geospatial_Mutation_Matrix.csv")

    print("=" * 60)
    print(f"Sanitized metadata genomes: {len(meta)}")
    print(f"Regulatory mutation rows: {len(genomics)}")
    print(f"Merged rows: {len(merged)}")
    print(f"Unique locations: {merged['geographic_location_clean'].nunique()}")
    print(f"Unique mutations: {merged['Mutation'].nunique()}")
    print(f"Output dir: {out_dir}")
    print("=" * 60)


def cmd_geospatial_heatmap(args):
    import matplotlib.pyplot as plt
    import seaborn as sns

    input_csv = Path(args.input_csv)
    out_matrix = Path(args.output_matrix_csv)
    out_plot = Path(args.output_plot)
    out_matrix.parent.mkdir(parents=True, exist_ok=True)
    out_plot.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_csv, dtype=str, keep_default_na=False)
    for c in ["Gene", "Mutation", "geographic_location_clean"]:
        if c not in df.columns:
            raise KeyError(f"Missing required column: {c}")

    df["Gene"] = df["Gene"].astype(str).str.strip()
    df["Mutation"] = df["Mutation"].astype(str).str.strip()
    df["Gene_Mutation"] = df["Gene"] + ":" + df["Mutation"]
    df["geographic_location_clean"] = df["geographic_location_clean"].astype(str).str.strip()
    df.loc[
        df["geographic_location_clean"].isin(["", "nan", "None"]),
        "geographic_location_clean",
    ] = pd.NA
    df = df.dropna(subset=["geographic_location_clean"]).copy()

    exclude = {x.strip().lower() for x in args.exclude_locations}
    df = df[~df["geographic_location_clean"].str.lower().isin(exclude)].copy()
    df = df[(df["Gene"] != "") & (df["Mutation"] != "") & (df["Gene_Mutation"] != ":")]

    freq = df.groupby(["geographic_location_clean", "Gene_Mutation"], dropna=False).size().reset_index(name="mutation_frequency")
    matrix = freq.pivot(index="geographic_location_clean", columns="Gene_Mutation", values="mutation_frequency").fillna(0).astype(int)
    matrix.to_csv(out_matrix)

    top_n = max(1, args.top_n)
    top_rows = matrix.sum(axis=1).sort_values(ascending=False).head(top_n).index
    top_cols = matrix.sum(axis=0).sort_values(ascending=False).head(top_n).index
    top = matrix.loc[top_rows, top_cols]

    plt.figure(figsize=(14, 10))
    ax = sns.heatmap(top, cmap="Reds", annot=True, fmt=".0f", linewidths=0.5, cbar_kws={"label": "Mutation Frequency"})
    ax.set_title("Geospatial Distribution of Regulatory AMR Mutations Across Top Clinical Sites", fontsize=14, pad=12)
    ax.set_xlabel("Gene:Mutation", fontsize=12)
    ax.set_ylabel("Clinical Location", fontsize=12)
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(out_plot, dpi=300)
    plt.close()

    print("=" * 60)
    print(f"Filtered rows: {len(df)}")
    print(f"Unique locations: {df['geographic_location_clean'].nunique()}")
    print(f"Unique Gene_Mutation labels: {df['Gene_Mutation'].nunique()}")
    print(f"Matrix saved: {out_matrix}")
    print(f"Heatmap saved: {out_plot}")
    print("=" * 60)


def _n50(lengths):
    if not lengths:
        return 0
    s = sorted(lengths, reverse=True)
    half = sum(s) / 2
    run = 0
    for x in s:
        run += x
        if run >= half:
            return x
    return 0


def cmd_qc_genomes(args):
    input_dir = Path(args.input_dir)
    out_csv = Path(args.output_summary_csv)
    ready_dir = Path(args.output_ready_dir)
    ready_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    contig_counter = Counter()
    files = sorted(input_dir.glob("*.fna"))

    for fp in files:
        contig_lengths = []
        total_bp = 0
        n_count = 0
        cur = 0
        with fp.open("r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith(">"):
                    if cur:
                        contig_lengths.append(cur)
                        cur = 0
                    continue
                seq = line.strip().upper()
                if not seq:
                    continue
                L = len(seq)
                cur += L
                total_bp += L
                n_count += seq.count("N")
        if cur:
            contig_lengths.append(cur)

        contigs = len(contig_lengths)
        contig_counter[contigs] += 1
        n50v = _n50(contig_lengths)
        n_frac = (n_count / total_bp) if total_bp else 1.0

        flags = []
        if total_bp < args.min_bp:
            flags.append("small_assembly")
        if total_bp > args.max_bp:
            flags.append("large_assembly")
        if contigs > args.high_fragmentation:
            flags.append("high_fragmentation")
        elif contigs > args.moderate_fragmentation:
            flags.append("moderate_fragmentation")
        if n50v < args.low_n50:
            flags.append("low_n50")
        elif n50v < args.moderate_n50:
            flags.append("moderate_n50")
        if n_frac > args.high_n_fraction:
            flags.append("high_ambiguous_bases")
        elif n_frac > args.moderate_n_fraction:
            flags.append("moderate_ambiguous_bases")

        severe = {"small_assembly", "high_fragmentation", "low_n50", "high_ambiguous_bases"}
        moderate = {"large_assembly", "moderate_fragmentation", "moderate_n50", "moderate_ambiguous_bases"}
        if any(f in severe for f in flags):
            quality = "poor"
        elif any(f in moderate for f in flags):
            quality = "moderate"
        else:
            quality = "good"

        ready = (
            args.min_bp <= total_bp <= args.max_bp
            and contigs <= args.ready_max_contigs
            and n50v >= args.ready_min_n50
            and n_frac <= args.ready_max_n_fraction
        )

        rows.append(
            {
                "genome_id": fp.stem,
                "file_name": fp.name,
                "total_bp": total_bp,
                "contigs": contigs,
                "n50": n50v,
                "n_fraction": round(n_frac, 6),
                "quality": quality,
                "extraction_ready": int(ready),
                "flags": ";".join(flags) if flags else "none",
            }
        )
        if ready:
            shutil.copy2(fp, ready_dir / fp.name)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ["genome_id"])
        w.writeheader()
        w.writerows(rows)

    q = Counter(r["quality"] for r in rows)
    print("=" * 60)
    print(f"Analyzed genomes: {len(rows)}")
    print(f"Quality counts: good={q.get('good',0)} moderate={q.get('moderate',0)} poor={q.get('poor',0)}")
    print(f"Extraction-ready genomes: {sum(r['extraction_ready'] for r in rows)}")
    print(f"Summary written: {out_csv}")
    print(f"Extraction-ready folder: {ready_dir}")
    print("=" * 60)


def cmd_presentation_plots(args):
    import matplotlib.pyplot as plt
    import networkx as nx
    import numpy as np
    import seaborn as sns
    from matplotlib.patches import Patch

    runs = args.runs
    base_output_dir = Path(args.base_output_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    regulatory_genes = {"marr", "acrr"}
    structural_genes = {"acra", "acrb", "tolc"}
    regulatory_color = "crimson"
    structural_color = "steelblue"
    other_color = "gray"

    sns.set_theme(style="whitegrid", context="talk")

    def safe_read_csv(csv_path):
        try:
            if not csv_path.exists() or csv_path.stat().st_size == 0:
                print(f"[WARN] Missing or empty file: {csv_path}")
                return pd.DataFrame()
            return pd.read_csv(csv_path)
        except Exception as exc:
            print(f"[WARN] Failed to read {csv_path}: {exc}")
            return pd.DataFrame()

    def parse_gene(node):
        token = str(node or "").strip()
        if ":" in token:
            return token.split(":", 1)[0].strip().lower()
        return ""

    def classify(gene):
        if gene in regulatory_genes:
            return "regulatory"
        if gene in structural_genes:
            return "structural"
        return "other"

    def color(cls):
        if cls == "regulatory":
            return regulatory_color
        if cls == "structural":
            return structural_color
        return other_color

    def scale(vals, lo, hi):
        if not vals:
            return []
        arr = np.asarray(vals, dtype=float)
        arr = np.nan_to_num(arr, nan=float(np.nanmin(arr)) if np.isfinite(np.nanmin(arr)) else 0.0)
        vmin = float(np.min(arr))
        vmax = float(np.max(arr))
        if vmax - vmin < 1e-12:
            return [0.5 * (lo + hi)] * len(arr)
        return (lo + (arr - vmin) * (hi - lo) / (vmax - vmin)).tolist()

    def build_graph(df):
        g = nx.Graph()
        node_freq = {}
        edge_w = {}
        if df.empty:
            return g, node_freq, edge_w
        required = {"Node_1", "Node_2", "Frequency", "Composite_Network_Score"}
        if not required.issubset(df.columns):
            print("[WARN] Epistasis CSV missing required columns; skipping graph content.")
            return g, node_freq, edge_w
        for _, row in df.iterrows():
            n1 = str(row.get("Node_1", "")).strip()
            n2 = str(row.get("Node_2", "")).strip()
            if not n1 or not n2:
                continue
            freq = float(row.get("Frequency", 0) or 0)
            score = float(row.get("Composite_Network_Score", 0) or 0)
            g.add_edge(n1, n2)
            node_freq[n1] = node_freq.get(n1, 0.0) + freq
            node_freq[n2] = node_freq.get(n2, 0.0) + freq
            ek = tuple(sorted((n1, n2)))
            edge_w[ek] = max(edge_w.get(ek, 0.0), score)
        return g, node_freq, edge_w

    def flatten_freq(epi_df):
        if epi_df.empty:
            return pd.DataFrame(columns=["mutation", "frequency", "gene", "gene_class"])
        needed = {"Node_1", "Node_2", "Frequency"}
        if not needed.issubset(epi_df.columns):
            return pd.DataFrame(columns=["mutation", "frequency", "gene", "gene_class"])
        rows = []
        for _, row in epi_df.iterrows():
            freq = float(row.get("Frequency", 0) or 0)
            for col in ("Node_1", "Node_2"):
                mut = str(row.get(col, "")).strip()
                if not mut:
                    continue
                gene = parse_gene(mut)
                rows.append({"mutation": mut, "frequency": freq, "gene": gene, "gene_class": classify(gene)})
        if not rows:
            return pd.DataFrame(columns=["mutation", "frequency", "gene", "gene_class"])
        out = pd.DataFrame(rows)
        out = out.groupby(["mutation", "gene", "gene_class"])["frequency"].sum().reset_index()
        return out.sort_values(by=["frequency"], ascending=[False])

    def load_biophysics(path, run_name):
        df = safe_read_csv(path)
        if df.empty:
            return pd.DataFrame(columns=["run", "mutation_network", "delta_delta_g", "status", "interpretation"])
        for c, default in [("mutation_network", None), ("delta_delta_g", np.nan), ("status", "unknown"), ("interpretation", None)]:
            if c not in df.columns:
                df[c] = default
        out = df[["mutation_network", "delta_delta_g", "status", "interpretation"]].copy()
        out["run"] = run_name
        out["delta_delta_g"] = pd.to_numeric(out["delta_delta_g"], errors="coerce")
        return out

    def triage(row):
        status = str(row.get("status", "")).strip()
        interp = str(row.get("interpretation", "")).strip()
        ddg = row.get("delta_delta_g", np.nan)
        if "FAILED_QC" in status or "failed" in status.lower():
            return "FAILED_QC_or_DockingFailed"
        if interp == "NO_SIGNIFICANT_BINDING_CHANGE":
            return "NO_SIGNIFICANT_BINDING_CHANGE"
        if pd.notna(ddg) and abs(float(ddg)) < 0.5:
            return "NearZero_DDG"
        return "Shift_or_Other"

    epistasis_by_run = {}
    biophysics_frames = []
    for run_name in runs:
        run_dir = base_output_dir / run_name
        _ = safe_read_csv(run_dir / "1_genomics_report.csv")
        epistasis_by_run[run_name] = safe_read_csv(run_dir / "2_epistasis_networks.csv")
        biophysics_frames.append(load_biophysics(run_dir / "3_biophysics_docking.csv", run_name))

    biophysics_all = pd.concat(biophysics_frames, ignore_index=True) if biophysics_frames else pd.DataFrame()

    # Plot 1: Epistasis web 3-panel.
    fig, axes = plt.subplots(1, len(runs), figsize=(8 * max(1, len(runs)), 8), constrained_layout=True)
    if len(runs) == 1:
        axes = [axes]
    legend_items = [
        Patch(facecolor=regulatory_color, edgecolor="none", label="Regulatory (marR, acrR)"),
        Patch(facecolor=structural_color, edgecolor="none", label="Structural (acrA, acrB, tolC)"),
        Patch(facecolor=other_color, edgecolor="none", label="Other"),
    ]
    for ax, run_name in zip(axes, runs):
        df = epistasis_by_run.get(run_name, pd.DataFrame())
        if df.empty:
            ax.axis("off")
            ax.text(0.5, 0.5, f"{run_name}\\nNo epistasis data", ha="center", va="center", fontsize=13)
            continue
        g, node_freq, edge_scores = build_graph(df)
        if g.number_of_nodes() == 0:
            ax.axis("off")
            ax.text(0.5, 0.5, f"{run_name}\\nNo valid network edges", ha="center", va="center", fontsize=13)
            continue
        pos = nx.spring_layout(g, seed=42, k=0.9 / max(g.number_of_nodes(), 2))
        nodes = list(g.nodes())
        edges = list(g.edges())
        node_sizes = scale([node_freq.get(n, 1.0) for n in nodes], 250, 2200)
        node_colors = [color(classify(parse_gene(n))) for n in nodes]
        edge_vals = [edge_scores.get(tuple(sorted(e)), 0.0) for e in edges]
        edge_widths = scale(edge_vals, 0.8, 5.0) if edge_vals else []
        nx.draw_networkx_edges(g, pos, ax=ax, width=edge_widths, alpha=0.45, edge_color="#444444")
        nx.draw_networkx_nodes(g, pos, ax=ax, node_size=node_sizes, node_color=node_colors, alpha=0.9)
        nx.draw_networkx_labels(g, pos, ax=ax, font_size=8)
        ax.set_title(run_name.replace("_", " "), fontsize=15, fontweight="bold")
        ax.set_axis_off()
    fig.suptitle("Epistasis Web Across Bulky Antibiotic Runs", fontsize=20, fontweight="bold")
    fig.legend(handles=legend_items, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.04))
    fig.savefig(out_dir / "epistasis_web_3panel.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Plot 2: Regulatory dominance top-15.
    fig, axes = plt.subplots(1, len(runs), figsize=(8 * max(1, len(runs)), 10), constrained_layout=True)
    if len(runs) == 1:
        axes = [axes]
    for ax, run_name in zip(axes, runs):
        flat = flatten_freq(epistasis_by_run.get(run_name, pd.DataFrame()))
        if flat.empty:
            ax.axis("off")
            ax.text(0.5, 0.5, f"{run_name}\\nNo mutation frequencies", ha="center", va="center", fontsize=13)
            continue
        top15 = flat.head(15).copy().sort_values("frequency", ascending=True)
        colors = [color(c) for c in top15["gene_class"]]
        ax.barh(top15["mutation"], top15["frequency"], color=colors, alpha=0.9)
        ax.set_title(run_name.replace("_", " "), fontsize=15, fontweight="bold")
        ax.set_xlabel("Aggregated Frequency")
        ax.set_ylabel("Mutation")
        ax.tick_params(axis="y", labelsize=9)
    fig.suptitle("Regulatory Dominance: Top 15 Mutation Frequencies", fontsize=20, fontweight="bold")
    fig.legend(handles=legend_items, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.03))
    fig.savefig(out_dir / "regulatory_dominance_top15.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Plot 3: MVBM biophysics triage summary.
    fig, axes = plt.subplots(1, 2, figsize=(18, 7), constrained_layout=True)
    if biophysics_all.empty:
        for ax in axes:
            ax.axis("off")
            ax.text(0.5, 0.5, "No biophysics data", ha="center", va="center", fontsize=14)
        fig.suptitle("MVBM Biophysics Triage Summary", fontsize=18, fontweight="bold")
        fig.savefig(out_dir / "mvbm_biophysics_triage_summary.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
    else:
        plot_df = biophysics_all.copy()
        plot_df["triage_category"] = plot_df.apply(triage, axis=1)
        count_df = plot_df.groupby(["run", "triage_category"], as_index=False).size().rename(columns={"size": "count"})
        sns.barplot(data=count_df, x="run", y="count", hue="triage_category", ax=axes[0])
        axes[0].set_title("Biophysics Outcome Categories")
        axes[0].set_xlabel("Run")
        axes[0].set_ylabel("Count")
        axes[0].tick_params(axis="x", rotation=15)
        sns.stripplot(data=plot_df, x="run", y="delta_delta_g", hue="triage_category", dodge=False, alpha=0.85, size=7, ax=axes[1])
        axes[1].axhline(0.0, color="black", linewidth=1.0, linestyle="--")
        axes[1].axhline(0.5, color="#555555", linewidth=0.9, linestyle=":")
        axes[1].axhline(-0.5, color="#555555", linewidth=0.9, linestyle=":")
        axes[1].set_title("$\\Delta\\Delta G$ Distribution by Run")
        axes[1].set_xlabel("Run")
        axes[1].set_ylabel("$\\Delta\\Delta G$ (kcal/mol)")
        axes[1].tick_params(axis="x", rotation=15)
        for ax in axes:
            handles, labels = ax.get_legend_handles_labels()
            if handles and ax.legend_ is not None:
                ax.legend_.remove()
        handles, labels = axes[0].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.02))
        fig.suptitle("MVBM Biophysics Triage Summary", fontsize=18, fontweight="bold")
        fig.savefig(out_dir / "mvbm_biophysics_triage_summary.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

        summary = (
            plot_df.groupby(["run", "triage_category"], as_index=False)
            .size()
            .rename(columns={"size": "count"})
            .sort_values(["run", "count"], ascending=[True, False])
        )
        summary.to_csv(out_dir / "mvbm_biophysics_triage_counts.csv", index=False)

    print(f"Presentation plots generated in: {out_dir}")


def build_parser():
    p = argparse.ArgumentParser(description="MutationScan supporting toolkit (consolidated utility operations).")
    sp = p.add_subparsers(dest="cmd", required=True)

    d = sp.add_parser("download-rest", help="Download BV-BRC genomes via REST API")
    d.add_argument("--csv-file", required=True)
    d.add_argument("--output-dir", required=True)
    d.add_argument("--threads", type=int, default=8)
    d.add_argument("--timeout", type=int, default=30)
    d.add_argument("--retries", type=int, default=3)
    d.add_argument("--min-bytes", type=int, default=1000)
    d.add_argument("--genome-id-column", default="Genome ID")
    d.add_argument("--failed-log", default="failed_rest_api_ids.txt")
    d.add_argument("--missing-log", default="missing_after_rest_api.txt")
    d.add_argument("--api-limit", type=int, default=5000)
    d.set_defaults(func=cmd_download_rest)

    m = sp.add_parser("fetch-metadata", help="Fetch BV-BRC genome metadata and enrich CSV")
    m.add_argument("--input-csv", required=True)
    m.add_argument("--genome-id-column", default="Genome ID")
    m.add_argument("--output-metadata-csv", required=True)
    m.add_argument("--output-enriched-csv", required=True)
    m.add_argument("--failed-log", required=True)
    m.add_argument("--batch-size", type=int, default=80)
    m.add_argument("--timeout", type=int, default=45)
    m.add_argument("--retries", type=int, default=3)
    m.add_argument("--sleep-seconds", type=float, default=0.1)
    m.set_defaults(func=cmd_fetch_metadata)

    g = sp.add_parser("geospatial-matrix", help="Build regulatory geospatial mutation matrix")
    g.add_argument("--metadata-csvs", nargs="+", required=True)
    g.add_argument("--genomics-reports", nargs="+", required=True)
    g.add_argument("--regulatory-genes", nargs="+", default=["marR", "acrR"])
    g.add_argument("--output-dir", required=True)
    g.set_defaults(func=cmd_geospatial_matrix)

    h = sp.add_parser("geospatial-heatmap", help="Generate top-N geospatial heatmap")
    h.add_argument("--input-csv", required=True)
    h.add_argument("--output-matrix-csv", required=True)
    h.add_argument("--output-plot", required=True)
    h.add_argument("--top-n", type=int, default=15)
    h.add_argument("--exclude-locations", nargs="+", default=["Asia", "Egypt", "United Kingdom", "UK", "Unknown", "UNKNOWN", "India"])
    h.set_defaults(func=cmd_geospatial_heatmap)

    q = sp.add_parser("qc-genomes", help="Run assembly QC and build extraction-ready subset")
    q.add_argument("--input-dir", required=True)
    q.add_argument("--output-summary-csv", required=True)
    q.add_argument("--output-ready-dir", required=True)
    q.add_argument("--min-bp", type=int, default=4_000_000)
    q.add_argument("--max-bp", type=int, default=6_500_000)
    q.add_argument("--moderate-fragmentation", type=int, default=200)
    q.add_argument("--high-fragmentation", type=int, default=500)
    q.add_argument("--low-n50", type=int, default=20_000)
    q.add_argument("--moderate-n50", type=int, default=50_000)
    q.add_argument("--moderate-n-fraction", type=float, default=0.001)
    q.add_argument("--high-n-fraction", type=float, default=0.01)
    q.add_argument("--ready-max-contigs", type=int, default=300)
    q.add_argument("--ready-min-n50", type=int, default=50_000)
    q.add_argument("--ready-max-n-fraction", type=float, default=0.01)
    q.set_defaults(func=cmd_qc_genomes)

    pz = sp.add_parser("presentation-plots", help="Generate comparative presentation plots")
    pz.add_argument(
        "--runs",
        nargs="+",
        default=["Tetracycline_Run", "Chloramphenicol_Run", "Tigecycline_Run"],
        help="Run folder names under base output directory",
    )
    pz.add_argument("--base-output-dir", default="data/output")
    pz.add_argument("--output-dir", default="data/presentation")
    pz.set_defaults(func=cmd_presentation_plots)

    return p


def main():
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
