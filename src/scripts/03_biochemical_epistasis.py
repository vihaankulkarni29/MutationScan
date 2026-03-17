import logging
import os
import sys
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd

# Adjust path to access local modules
sys.path.append(str(Path(__file__).resolve().parents[2]))

try:
    from mutation_scan.analysis.control_scan import MutationScorer
except ImportError:
    # Fallback if installed as an external package
    from controlscan.scorer import MutationScorer

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Snakemake Context Injection
mutations_csv = Path(getattr(snakemake.input, "mutations_csv", snakemake.input.report))
epistasis_csv = Path(getattr(snakemake.output, "epistasis_csv", snakemake.output.networks))
networks_dir = Path(getattr(snakemake.output, "networks_dir", snakemake.output.plots_dir))
RESULTS_DIR = Path(snakemake.params.out_dir)
os.makedirs(RESULTS_DIR, exist_ok=True)

networks_dir.mkdir(parents=True, exist_ok=True)

logger.info("Step 2: Loading genomic report for biochemical scoring...")
df = pd.read_csv(mutations_csv)

if df.empty or 'Mutation' not in df.columns:
    logger.warning("No mutations to analyze. Generating empty Phase 2/3 outputs.")
    pd.DataFrame().to_csv(epistasis_csv, index=False)
    sys.exit(0)

# ---------------------------------------------------------
# PHASE 2: ControlScan Biochemical Scoring
# ---------------------------------------------------------
logger.info("Step 2.1: Applying ControlScan mathematical biochemist scoring...")
scorer = MutationScorer()


def apply_controlscan(mutation_str):
    try:
        # Assumes mutation format like 'S83L'
        result = scorer.score_single(mutation_str)
        severity = result.get('Severity')
        return float(severity) if severity is not None else 1.0
    except Exception:
        return 1.0  # Baseline fallback if parsing fails


df['ControlScan_Score'] = df['Mutation'].apply(apply_controlscan)

# ---------------------------------------------------------
# PHASE 3: Epistasis Networks & Composite Scoring
# ---------------------------------------------------------
logger.info("Step 3.1: Calculating Co-occurrence Frequencies...")
co_occurrences = {}
# Build a lookup dictionary for scores to speed up the loop
score_lookup = df.set_index('Mutation')['ControlScan_Score'].to_dict()

# Group by genome to find mutations that happen together
for accession, group in df.groupby('Accession'):
    muts = sorted(group['Mutation'].dropna().unique())
    if len(muts) > 1:
        for pair in combinations(muts, 2):
            if pair not in co_occurrences:
                co_occurrences[pair] = 0
            co_occurrences[pair] += 1

logger.info("Step 3.2: Applying Composite Epistatic Math...")
network_data = []
for pair, freq in co_occurrences.items():
    score_A = score_lookup.get(pair[0], 1.0)
    score_B = score_lookup.get(pair[1], 1.0)
    avg_severity = (score_A + score_B) / 2.0

    # The Core Equation
    composite_score = freq * avg_severity

    network_data.append(
        {
            'Node_1': pair[0],
            'Node_2': pair[1],
            'Frequency': freq,
            'Avg_Biochemical_Severity': avg_severity,
            'Composite_Network_Score': composite_score,
        }
    )

net_df = pd.DataFrame(network_data)
if not net_df.empty:
    net_df = net_df.sort_values('Composite_Network_Score', ascending=False).head(5)
net_df.to_csv(epistasis_csv, index=False)

# ---------------------------------------------------------
# PLOTTING: Generating the Visual Artifacts
# ---------------------------------------------------------
logger.info("Step 3.3: Generating Matplotlib Network Graphs...")
for rank, (_, row) in enumerate(net_df.iterrows(), start=1):
    G = nx.Graph()
    G.add_edge(row['Node_1'], row['Node_2'], weight=row['Composite_Network_Score'])

    plt.figure(figsize=(8, 6))
    pos = nx.spring_layout(G, seed=42)

    # Draw nodes and edges
    nx.draw_networkx_nodes(G, pos, node_color='#4A90E2', node_size=2500, edgecolors='black')
    nx.draw_networkx_edges(G, pos, width=3.0, edge_color='#505050')
    nx.draw_networkx_labels(G, pos, font_size=12, font_weight='bold', font_color='white')

    # Add Edge labels (The Score)
    edge_labels = {
        (row['Node_1'], row['Node_2']): f"Score: {row['Composite_Network_Score']:.1f}"
    }
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=10, font_color='red')

    plt.title(f"Rank {rank} Epistatic Network", fontsize=16, pad=20)
    plt.axis('off')

    plot_path = networks_dir / f"Network_Rank_{rank}.png"
    plt.savefig(plot_path, bbox_inches='tight', dpi=300)
    plt.close()

logger.info(f"Phase 2/3 Complete. Artifacts saved to {networks_dir}")
