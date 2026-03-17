import logging
import os
import re
import stat
import subprocess
import sys
import urllib.request
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

AA_MAP = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
}

SMINA_URL = "https://sourceforge.net/projects/smina/files/smina.static/download"
CIPRO_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/2764/SDF"
BOX_SIZE = 25.0
EXHAUSTIVENESS = 8


def parse_chain_map(raw):
    chain_map = {}
    for token in str(raw or "").split(","):
        token = token.strip()
        if not token or ":" not in token:
            continue
        gene, chain = token.split(":", 1)
        gene = gene.strip().lower()
        chain = chain.strip()
        if gene and chain:
            chain_map[gene] = chain
    return chain_map


def parse_mutation_token(token, chain_map, default_chain="A"):
    if not isinstance(token, str) or not token.strip():
        return None

    raw = token.strip()
    gene = None
    mutation = raw
    if ":" in raw:
        gene, mutation = raw.split(":", 1)
        gene = gene.strip().lower()

    match = re.match(r"^([A-Za-z])(\d+)([A-Za-z])$", mutation.strip())
    if not match:
        return None

    old_aa, residue_str, new_aa = match.groups()
    old_three = AA_MAP.get(old_aa.upper())
    new_three = AA_MAP.get(new_aa.upper())
    if not old_three or not new_three:
        return None

    residue = int(residue_str)
    chain = chain_map.get(gene, default_chain)
    pdbfixer_mut = f"{old_three}-{residue}-{new_three}"

    return {
        "gene": gene,
        "token": raw,
        "chain": chain,
        "residue": residue,
        "pdbfixer_mut": pdbfixer_mut,
    }


def find_ca_coord(pdb_path, chain_id, residue_num):
    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if not line.startswith("ATOM"):
                    continue
                if line[12:16].strip() != "CA":
                    continue
                if line[21].strip() != str(chain_id).strip():
                    continue
                try:
                    seq_num = int(line[22:26].strip())
                except ValueError:
                    continue
                if seq_num == int(residue_num):
                    return (
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    )
    except OSError:
        return None
    return None


def find_residue_name(pdb_path, chain_id, residue_num):
    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if not line.startswith("ATOM"):
                    continue
                if line[21].strip() != str(chain_id).strip():
                    continue
                try:
                    seq_num = int(line[22:26].strip())
                except ValueError:
                    continue
                if seq_num == int(residue_num):
                    return line[17:20].strip().upper()
    except OSError:
        return None
    return None


def pocket_center_from_mutations(pdb_path, parsed_mutations):
    coords = []
    for mut in parsed_mutations:
        coord = find_ca_coord(pdb_path, mut["chain"], mut["residue"])
        if coord is not None:
            coords.append(coord)

    if coords:
        n = len(coords)
        return (
            sum(c[0] for c in coords) / n,
            sum(c[1] for c in coords) / n,
            sum(c[2] for c in coords) / n,
        )

    fallback = find_ca_coord(pdb_path, "A", 596)
    if fallback is not None:
        return fallback
    return (18.5, 55.2, 20.1)


def run_cmd(cmd, label, allow_failure=False):
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0 and not allow_failure:
        logger.error("Command failed during %s", label)
        if result.stdout:
            logger.error(result.stdout)
        if result.stderr:
            logger.error(result.stderr)
    return result


def parse_affinity(text):
    for line in text.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0] == "1":
            try:
                return float(parts[1])
            except ValueError:
                continue
    return None


def ensure_smina_binary(work_dir):
    work_dir.mkdir(parents=True, exist_ok=True)
    smina_path = work_dir / "smina"
    if not smina_path.exists():
        logger.info("Downloading Smina binary to %s", smina_path)
        urllib.request.urlretrieve(SMINA_URL, smina_path)

    try:
        mode = os.stat(smina_path).st_mode
        os.chmod(smina_path, mode | stat.S_IEXEC)
    except OSError:
        pass

    qvina_candidate = Path(sys.executable).resolve().parent / "Library" / "bin" / "qvina2.exe"
    if os.name == "nt" and qvina_candidate.exists():
        logger.info("Using Windows qvina2 fallback at %s", qvina_candidate)
        return str(qvina_candidate)

    return str(smina_path)


def ensure_ligand_file(ligand_path):
    ligand_path.parent.mkdir(parents=True, exist_ok=True)
    if ligand_path.exists() and ligand_path.stat().st_size > 0:
        return
    logger.info("Ligand file missing, downloading %s", ligand_path)
    urllib.request.urlretrieve(CIPRO_URL, ligand_path)


def prepare_ligand_pdbqt(ligand_path, out_dir):
    ligand_pdbqt = out_dir / f"{ligand_path.stem}.pdbqt"
    if ligand_pdbqt.exists() and ligand_pdbqt.stat().st_size > 0:
        return ligand_pdbqt

    cmd = [
        sys.executable,
        "-m",
        "meeko.cli.mk_prepare_ligand",
        "-i",
        str(ligand_path),
        "-o",
        str(ligand_pdbqt),
    ]
    result = run_cmd(cmd, "ligand preparation", allow_failure=True)
    if result.returncode != 0 and not ligand_pdbqt.exists():
        raise RuntimeError("Ligand preparation failed and no PDBQT was generated.")
    return ligand_pdbqt


def prepare_receptor_pdbqt(receptor_pdb, center, out_prefix):
    receptor_pdbqt = Path(f"{out_prefix}.pdbqt")
    if receptor_pdbqt.exists() and receptor_pdbqt.stat().st_size > 0:
        return receptor_pdbqt

    base_cmd = [
        sys.executable,
        "-m",
        "meeko.cli.mk_prepare_receptor",
        "--read_pdb",
        str(receptor_pdb),
        "-o",
        str(out_prefix),
        "-p",
        "--box_center",
        str(center[0]),
        str(center[1]),
        str(center[2]),
        "--box_size",
        str(BOX_SIZE),
        str(BOX_SIZE),
        str(BOX_SIZE),
    ]

    strict = run_cmd(base_cmd, f"receptor preparation {receptor_pdb.name}", allow_failure=True)
    if strict.returncode == 0 and receptor_pdbqt.exists():
        return receptor_pdbqt

    logger.warning("Strict receptor preparation failed for %s. Retrying with allow_bad_res.", receptor_pdb)
    relaxed = run_cmd(base_cmd + ["-a"], f"receptor fallback {receptor_pdb.name}", allow_failure=True)
    if relaxed.returncode != 0 and not receptor_pdbqt.exists():
        raise RuntimeError(f"Receptor preparation failed for {receptor_pdb}")
    return receptor_pdbqt


def apply_mutations_to_structure(reference_pdb, parsed_mutations, mutated_pdb_out):
    import pdbfixer
    from openmm.app import PDBFile

    fixer = pdbfixer.PDBFixer(filename=str(reference_pdb))

    grouped = {}
    for mut in parsed_mutations:
        actual_old = find_residue_name(reference_pdb, mut["chain"], mut["residue"])
        if not actual_old:
            logger.warning(
                "Could not resolve residue identity for %s at %s:%s. Skipping this mutation token.",
                mut["token"],
                mut["chain"],
                mut["residue"],
            )
            continue
        requested_new = mut["pdbfixer_mut"].split("-")[-1]
        reconciled = f"{actual_old}-{mut['residue']}-{requested_new}"
        grouped.setdefault(mut["chain"], []).append(reconciled)

    if not grouped:
        raise RuntimeError("No valid mutations could be reconciled against reference PDB residues.")

    for chain_id, mut_list in grouped.items():
        fixer.applyMutations(mut_list, chain_id)

    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    with open(mutated_pdb_out, "w", encoding="utf-8") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle)


def run_docking(binary, receptor_pdbqt, ligand_pdbqt, center, out_pose):
    cmd = [
        binary,
        "--receptor",
        str(receptor_pdbqt),
        "--ligand",
        str(ligand_pdbqt),
        "--center_x",
        str(center[0]),
        "--center_y",
        str(center[1]),
        "--center_z",
        str(center[2]),
        "--size_x",
        str(BOX_SIZE),
        "--size_y",
        str(BOX_SIZE),
        "--size_z",
        str(BOX_SIZE),
        "--exhaustiveness",
        str(EXHAUSTIVENESS),
        "--out",
        str(out_pose),
    ]
    result = run_cmd(cmd, f"docking {receptor_pdbqt.name}", allow_failure=True)
    if result.returncode != 0:
        return None
    return parse_affinity(result.stdout)


def append_disclaimer(readme_path, protein_name, ligand_name, mutation_network, wt_affinity, mut_affinity, ddg_score):
    block = f"""=====================================================================
MUTATIONSCAN: BEST-EFFORT PARTIAL DOCKING SIMULATION REPORT
=====================================================================
TARGET: {protein_name}
LIGAND: {ligand_name}
EPISTATIC NETWORK: {mutation_network}

RESULTS:
- Wild-Type Baseline Affinity: {wt_affinity} kcal/mol
- Mutated Protein Affinity: {mut_affinity} kcal/mol
- Thermodynamic Shift (ΔΔG): {ddg_score} kcal/mol

METHODOLOGICAL DISCLAIMER:
What we have conducted is a targeted, high-throughput partial docking simulation for the mutation(s) listed above. This biophysical module utilizes a rigid-receptor, flexible-ligand local topology search via the Smina (AutoDock Vina) scoring function. 

While this rapid geometric search is highly effective for estimating relative thermodynamic shifts (ΔΔG) and localized steric clashes, it is a "best-effort" partial docking. It does not account for large-scale backbone conformational flexibility, explicit solvent dynamics, or microsecond-scale entropic effects. For absolute free-energy perturbation (FEP) calculations or fully accurate binding affinities, please conduct a full-scale Molecular Dynamics (MD) simulation using a suite such as Gromacs or AMBER coupled with thermodynamic integration.
=====================================================================
"""
    with open(readme_path, "a", encoding="utf-8") as handle:
        handle.write(block)


def fmt_number(value):
    if value is None:
        return "N/A"
    return f"{value:.3f}"


def main():
    RESULTS_DIR = Path(snakemake.params.out_dir)
    os.makedirs(RESULTS_DIR, exist_ok=True)

    networks_csv = Path(snakemake.input.networks)
    reference_pdb = Path(snakemake.input.pdb_file)
    docking_report = Path(snakemake.output.docking_report)
    mutated_pdbs_dir = Path(snakemake.output.mutated_pdbs)
    readme_path = Path(snakemake.output[2])

    ligand_cfg = snakemake.config.get("ligand", "")
    ligand_path = Path(ligand_cfg) if ligand_cfg else RESULTS_DIR / "ciprofloxacin.sdf"
    chain_map = parse_chain_map(snakemake.config.get("chain_map", ""))
    default_chain = "A"

    mutated_pdbs_dir.mkdir(parents=True, exist_ok=True)
    docking_report.parent.mkdir(parents=True, exist_ok=True)
    readme_path.parent.mkdir(parents=True, exist_ok=True)
    readme_path.write_text("", encoding="utf-8")

    if not networks_csv.exists():
        raise FileNotFoundError(f"Missing epistasis network file: {networks_csv}")
    if not reference_pdb.exists():
        raise FileNotFoundError(f"Missing reference PDB file: {reference_pdb}")

    ensure_ligand_file(ligand_path)
    docking_binary = ensure_smina_binary(mutated_pdbs_dir)
    ligand_pdbqt = prepare_ligand_pdbqt(ligand_path, mutated_pdbs_dir)

    df = pd.read_csv(networks_csv)
    if df.empty or not {"Node_1", "Node_2"}.issubset(df.columns):
        pd.DataFrame(columns=[
            "Node_1",
            "Node_2",
            "mutation_network",
            "wt_affinity",
            "mut_affinity",
            "delta_delta_g",
            "status",
            "center_x",
            "center_y",
            "center_z",
        ]).to_csv(docking_report, index=False)
        readme_path.write_text("No valid epistatic network rows found.\n", encoding="utf-8")
        return

    results = []
    for idx, row in df.iterrows():
        node_1 = str(row["Node_1"]).strip()
        node_2 = str(row["Node_2"]).strip()
        mutation_network = f"{node_1} + {node_2}"

        parsed_mutations = []
        for token in (node_1, node_2):
            parsed = parse_mutation_token(token, chain_map, default_chain=default_chain)
            if parsed is not None:
                parsed_mutations.append(parsed)

        if not parsed_mutations:
            results.append(
                {
                    "Node_1": node_1,
                    "Node_2": node_2,
                    "mutation_network": mutation_network,
                    "wt_affinity": None,
                    "mut_affinity": None,
                    "delta_delta_g": None,
                    "status": "skipped_unparseable_mutations",
                    "center_x": None,
                    "center_y": None,
                    "center_z": None,
                }
            )
            append_disclaimer(
                readme_path,
                reference_pdb.stem,
                ligand_path.name,
                mutation_network,
                "N/A",
                "N/A",
                "N/A",
            )
            continue

        center = pocket_center_from_mutations(reference_pdb, parsed_mutations)

        wt_receptor_prefix = mutated_pdbs_dir / f"network_{idx + 1}_WT_receptor"
        wt_receptor_pdbqt = prepare_receptor_pdbqt(reference_pdb, center, wt_receptor_prefix)
        wt_pose = mutated_pdbs_dir / f"network_{idx + 1}_WT_docked.pdbqt"
        wt_affinity = run_docking(docking_binary, wt_receptor_pdbqt, ligand_pdbqt, center, wt_pose)

        mutated_pdb = mutated_pdbs_dir / f"network_{idx + 1}_mutated.pdb"
        mut_affinity = None
        status = "ok"
        try:
            apply_mutations_to_structure(reference_pdb, parsed_mutations, mutated_pdb)
            mut_receptor_prefix = mutated_pdbs_dir / f"network_{idx + 1}_MUT_receptor"
            mut_receptor_pdbqt = prepare_receptor_pdbqt(mutated_pdb, center, mut_receptor_prefix)
            mut_pose = mutated_pdbs_dir / f"network_{idx + 1}_MUT_docked.pdbqt"
            mut_affinity = run_docking(docking_binary, mut_receptor_pdbqt, ligand_pdbqt, center, mut_pose)
        except Exception as exc:
            logger.warning("Mutation threading or mutant docking failed for %s: %s", mutation_network, exc)
            status = "mutant_failed_fallback"

        ddg = None
        if wt_affinity is not None and mut_affinity is not None:
            ddg = mut_affinity - wt_affinity

        if wt_affinity is None:
            status = "wildtype_docking_failed"

        results.append(
            {
                "Node_1": node_1,
                "Node_2": node_2,
                "mutation_network": mutation_network,
                "wt_affinity": wt_affinity,
                "mut_affinity": mut_affinity,
                "delta_delta_g": ddg,
                "status": status,
                "center_x": center[0],
                "center_y": center[1],
                "center_z": center[2],
            }
        )

        append_disclaimer(
            readme_path,
            reference_pdb.stem,
            ligand_path.name,
            mutation_network,
            fmt_number(wt_affinity),
            fmt_number(mut_affinity),
            fmt_number(ddg),
        )

    pd.DataFrame(results).to_csv(docking_report, index=False)


if __name__ == "__main__":
    main()
