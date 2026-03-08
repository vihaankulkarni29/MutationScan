import logging
import re
import sys
from pathlib import Path

import pandas as pd

# Adjust path to ensure the mutation_scan package is discoverable
sys.path.append(str(Path(__file__).resolve().parents[2]))

try:
    from mutation_scan.biophysics.autoscan_bridge import AutoScanBridge
except ImportError:
    logging.error(
        "CRITICAL: AutoScanBridge module not found or missing biophysics dependencies."
    )
    sys.exit(1)

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# Snakemake Context Injection
epistasis_csv = Path(snakemake.input.epistasis_csv)
docking_csv = Path(snakemake.output.docking_csv)
mutated_dir = Path(snakemake.output.mutated_dir)
readme_path = Path(snakemake.output.readme)

pdb_file = str(snakemake.params.pdb)
chain_map_str = str(snakemake.params.chain_map)
ligand_file = str(snakemake.params.ligand)
cx = float(snakemake.params.center_x)
cy = float(snakemake.params.center_y)
cz = float(snakemake.params.center_z)
stiffness = float(snakemake.params.stiffness)

mutated_dir.mkdir(parents=True, exist_ok=True)


def parse_chain_map(raw: str) -> dict:
    chain_map = {}
    if not raw:
        return chain_map

    for mapping in raw.split(","):
        mapping = mapping.strip()
        if not mapping:
            continue
        if ":" not in mapping:
            logger.warning("Skipping invalid chain_map token: %s", mapping)
            continue
        gene, chain = mapping.split(":", 1)
        chain_map[gene.strip().lower()] = chain.strip()
    return chain_map


def extract_residue_number(mutation: str):
    if not isinstance(mutation, str):
        return None
    match = re.match(r"^[A-Za-z](\d+)[A-Za-z*]$", mutation.strip())
    if not match:
        return None
    return int(match.group(1))


def parse_gene_mutation(token: str):
    """Parse gene-prefixed mutation tokens like 'gyrA:S2L'."""
    if not isinstance(token, str):
        return None, None

    token = token.strip()
    if not token:
        return None, None

    if ":" in token:
        gene, mutation = token.split(":", 1)
        residue = extract_residue_number(mutation)
        return gene.strip().lower(), residue

    return None, extract_residue_number(token)


def load_residue_index_from_pdb(pdb_path: str, chain: str) -> set:
    """Return residue numbers that have ATOM coordinates for the requested chain."""
    residues = set()
    pdb_file_path = Path(pdb_path)
    if not pdb_file_path.exists():
        return residues

    try:
        with pdb_file_path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if not line.startswith("ATOM"):
                    continue
                line_chain = line[21].strip() if len(line) > 21 else ""
                if line_chain != chain:
                    continue
                seq_id_raw = line[22:26].strip() if len(line) >= 26 else ""
                if not seq_id_raw:
                    continue
                try:
                    residues.add(int(seq_id_raw))
                except ValueError:
                    continue
    except OSError as exc:
        logger.warning("Could not read PDB for residue pre-check: %s", exc)

    return residues


def expected_chain_for_gene(gene: str, default_chain_value: str) -> str:
    if gene and gene in chain_map:
        return chain_map[gene]
    return default_chain_value


chain_map = parse_chain_map(chain_map_str)
default_chain = next(iter(chain_map.values()), "A")

if not pdb_file:
    logger.error("Phase 4 requested without default_pdb. Set config default_pdb to enable HTVS.")
    sys.exit(1)

if not epistasis_csv.exists():
    logger.error("Missing epistasis input: %s", epistasis_csv)
    sys.exit(1)

logger.info("Step 4.0: Loading epistasis networks from %s", epistasis_csv)
try:
    df = pd.read_csv(epistasis_csv)
except pd.errors.EmptyDataError:
    logger.warning("Epistasis CSV is empty (no columns). Generating empty Phase 4 outputs.")
    pd.DataFrame().to_csv(docking_csv, index=False)
    readme_path.write_text("No epistatic network data to process.\n", encoding="utf-8")
    sys.exit(0)

required_columns = {"Node_1", "Node_2"}
if df.empty or not required_columns.issubset(df.columns):
    logger.warning("No epistatic networks found. Generating empty Phase 4 outputs.")
    pd.DataFrame().to_csv(docking_csv, index=False)
    readme_path.write_text("No epistatic network data to process.\n", encoding="utf-8")
    sys.exit(0)

networks = df[["Node_1", "Node_2"]].dropna().head(5).values.tolist()

logger.info("Step 4.1: Initializing AutoScanBridge with baseline PDB: %s", pdb_file)
try:
    # Newer bridge APIs may accept explicit docking/minimization parameters.
    bridge = AutoScanBridge(
        default_pdb=pdb_file,
        chain_map=chain_map,
        ligand_sdf=ligand_file,
        output_dir=mutated_dir,
        center=(cx, cy, cz),
        stiffness=stiffness,
    )
except TypeError:
    # Backward-compatible fallback to empirical proxy API currently in repo.
    bridge = AutoScanBridge()
    logger.info("AutoScanBridge running in compatibility mode (empirical proxy backend).")

logger.info("Step 4.2: Running HTVS biophysics pipeline over top %d networks...", len(networks))

try:
    if hasattr(bridge, "run_docking_pipeline"):
        results_df = bridge.run_docking_pipeline(networks)
        if not isinstance(results_df, pd.DataFrame):
            results_df = pd.DataFrame(results_df)
    else:
        rows = []
        pdb_id = Path(pdb_file).stem
        for node_1, node_2 in networks:
            gene_1, residue_1 = parse_gene_mutation(node_1)
            gene_2, residue_2 = parse_gene_mutation(node_2)

            chain_1 = expected_chain_for_gene(gene_1, default_chain)
            chain_2 = expected_chain_for_gene(gene_2, default_chain)
            chain = chain_1

            residues = [r for r in [residue_1, residue_2] if r is not None]
            residues = sorted(set(residues))

            if not residues:
                rows.append(
                    {
                        "Node_1": node_1,
                        "Node_2": node_2,
                        "PDB": pdb_id,
                        "Chain": chain,
                        "Residues": "",
                        "WT_Affinity": None,
                        "Mutant_Affinity": None,
                        "Delta_Delta_G": None,
                        "Status": "skipped_no_residue_parse",
                    }
                )
                continue

            missing_details = []
            if residue_1 is not None:
                available_chain_1 = load_residue_index_from_pdb(pdb_file, chain_1)
                if not available_chain_1 or residue_1 not in available_chain_1:
                    missing_details.append(f"{node_1} (chain {chain_1})")
            if residue_2 is not None:
                available_chain_2 = load_residue_index_from_pdb(pdb_file, chain_2)
                if not available_chain_2 or residue_2 not in available_chain_2:
                    missing_details.append(f"{node_2} (chain {chain_2})")

            if missing_details:
                logger.warning(
                    "AutoScanBridge pre-check: missing 3D coordinates for %s in %s. Skipping network.",
                    missing_details,
                    pdb_file,
                )
                rows.append(
                    {
                        "Node_1": node_1,
                        "Node_2": node_2,
                        "PDB": pdb_id,
                        "Chain": chain,
                        "Residues": ",".join(str(r) for r in residues),
                        "WT_Affinity": None,
                        "Mutant_Affinity": None,
                        "Delta_Delta_G": None,
                        "Status": "skipped_missing_3d_coordinates",
                    }
                )
                continue

            result = bridge.run_comparative_docking(
                pdb_id=pdb_id,
                chain=chain,
                residues=residues,
                ligand_smiles=ligand_file,
                output_dir=mutated_dir,
            )

            if not result:
                rows.append(
                    {
                        "Node_1": node_1,
                        "Node_2": node_2,
                        "PDB": pdb_id,
                        "Chain": chain,
                        "Residues": ",".join(str(r) for r in residues),
                        "WT_Affinity": None,
                        "Mutant_Affinity": None,
                        "Delta_Delta_G": None,
                        "Status": "failed",
                    }
                )
                continue

            rows.append(
                {
                    "Node_1": node_1,
                    "Node_2": node_2,
                    "PDB": pdb_id,
                    "Chain": default_chain,
                    "Residues": ",".join(str(r) for r in residues),
                    "WT_Affinity": result.get("wt_affinity"),
                    "Mutant_Affinity": result.get("mutant_affinity"),
                    "Delta_Delta_G": result.get("delta_delta_g"),
                    "Status": "ok",
                }
            )

        results_df = pd.DataFrame(rows)

    results_df.to_csv(docking_csv, index=False)
    logger.info("Phase 4 complete. Docking results saved to %s", docking_csv)
except Exception as e:
    logger.error("Biophysics execution failed: %s", e)
    sys.exit(1)

readme_content = """# MutationScan: HTVS Biophysics Report

Generated artifacts in this directory come from MutationScan Phase 4 HTVS biophysics.

This stage uses a fast triage path for mutated epistatic networks and produces DDG-style
comparative affinity metrics for downstream prioritization.

Call to action:
Treat these outputs as a ranking layer. For definitive structural thermodynamics, run full
molecular dynamics and free-energy workflows (for example GROMACS/AMBER + MMPBSA/FEP).
"""
readme_path.write_text(readme_content, encoding="utf-8")
