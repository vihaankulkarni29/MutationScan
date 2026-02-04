# MutationScan Configuration Directory

This directory contains organism-agnostic configuration files for customizing MutationScan behavior.

## Files

### drug_mapping.json
**Purpose:** Maps gene names to their associated antimicrobial drugs

**Format:**
```json
{
  "geneName": "drugName",
  "gyrA": "ciprofloxacin",
  "rpoB": "rifampicin"
}
```

**Usage:**
- Loaded by `VariantCaller` during initialization
- Determines resistance phenotype for detected mutations
- If a gene is not in this file, phenotype defaults to "Unknown"

**Customization:**
1. Add new gene-drug pairs for your organism
2. Update existing mappings for organism-specific nomenclature
3. Use lowercase gene names (case-insensitive matching applied)

---

### target_genes.txt
**Purpose:** Optional list of specific genes to analyze

**Format:**
```
# Comments start with #
geneName1
geneName2
acrA
acrB
```

**Usage:**
- Passed to `main.py` via `--targets` argument
- If not provided, GeneFinder uses default ABRicate broad screening
- Useful for focused studies (e.g., only efflux pump genes)

**Customization:**
1. List one gene name per line
2. Remove all lines to use ABRicate default (all resistance genes)
3. Use organism-specific gene names

---

## Usage Examples

### Example 1: E. coli Efflux Pump Study
```bash
# Create target_genes.txt with:
acrA
acrB
acrR
marR

# Run MutationScan
python src/main.py \
  --email your@email.com \
  --query "Escherichia coli" \
  --targets data/config/target_genes.txt \
  --config-dir data/config
```

### Example 2: M. tuberculosis First-Line Resistance
```bash
# Create target_genes.txt with:
rpoB
katG
inhA
embB
pncA

# Update drug_mapping.json with TB-specific drugs
{
  "rpoB": "rifampicin",
  "katG": "isoniazid",
  "inhA": "isoniazid",
  "embB": "ethambutol",
  "pncA": "pyrazinamide"
}

# Run MutationScan
python src/main.py \
  --email your@email.com \
  --query "Mycobacterium tuberculosis" \
  --targets data/config/target_genes.txt
```

### Example 3: Custom Organism (no hardcoded assumptions)
```bash
# 1. Create your own drug_mapping.json
{
  "customGeneA": "customDrugX",
  "customGeneB": "customDrugY"
}

# 2. Create target_genes.txt
customGeneA
customGeneB

# 3. Provide reference proteins in data/refs/
#    - customGeneA_WT.faa
#    - customGeneB_WT.faa

# 4. Run pipeline
python src/main.py \
  --email your@email.com \
  --query "Your Organism" \
  --targets data/config/target_genes.txt \
  --config-dir data/config
```

---

## Anti-Hardcoding Rules

**MutationScan is organism-agnostic.** All biological assumptions must be in configuration files:

✅ **Allowed in Python code:**
- Generic variable names: `gene`, `protein`, `mutation`
- Data structure definitions
- Algorithm logic
- File I/O operations

❌ **NOT allowed in Python code:**
- Specific gene names: `"gyrA"`, `"rpoB"`, `"acrR"`
- Organism names: `"Escherichia coli"`, `"Mycobacterium tuberculosis"`
- Drug names: `"ciprofloxacin"`, `"rifampicin"` (except in argparse defaults)
- Hardcoded reference sequences

**Exception:** Default values in argparse help strings are OK for documentation:
```python
parser.add_argument('--query', default='Escherichia coli',
                    help='NCBI search query (default: Escherichia coli)')
```

---

## Configuration Validation

MutationScan validates configuration files at runtime:

1. **drug_mapping.json**
   - Must be valid JSON
   - If missing, defaults to empty dict (logs warning)
   - Gene names are case-insensitive

2. **target_genes.txt**
   - Lines starting with `#` are comments
   - Empty lines ignored
   - If missing, GeneFinder uses ABRicate default

3. **ML Models (models/ directory)**
   - Automatically discovers `*.pkl` files
   - Model filename determines antibiotic (e.g., `ciprofloxacin_predictor.pkl`)
   - If no model matches, ML prediction disabled for that drug

---

## Extending MutationScan

To add support for a new organism or gene:

1. **Add gene-drug mapping** in `drug_mapping.json`
2. **Create reference protein** in `data/refs/{GeneName}_WT.faa`
3. **Optionally add to target list** in `target_genes.txt`
4. **Run pipeline** - no code changes needed!

**Example: Adding *blaNDM-1* gene**
```bash
# 1. Update drug_mapping.json
{
  "blaNDM-1": "carbapenems"
}

# 2. Add reference protein
# data/refs/blaNDM-1_WT.faa
>blaNDM-1_WT
MATALMTKLPF...

# 3. Add to target_genes.txt
blaNDM-1

# 4. Run pipeline (no code changes!)
python src/main.py --email user@example.com --query "Klebsiella pneumoniae"
```

Done! MutationScan automatically adapts to your organism.
