#!/usr/bin/env python3
"""
Demo Pipeline - Steps 2-5 for E. coli K-12 MG1655
(Step 1 already completed with direct FTP download)
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, '/app/src')

from mutation_scan.core import GeneFinder, SequenceExtractor
from mutation_scan.analysis import VariantCaller
from mutation_scan.visualization import PyMOLVisualizer

print("\n" + "="*70)
print("  MUTATIONSCAN DEMO - E. COLI K-12 EFFLUX PUMP RESISTANCE")
print("="*70)

# Paths
genomes_dir = Path("/app/data/genomes")
proteins_dir = Path("/app/data/proteins")
refs_dir = Path("/app/data/refs")
results_dir = Path("/app/data/results")
viz_dir = results_dir / "visualizations"

# Create directories
for d in [proteins_dir, refs_dir, results_dir, viz_dir]:
    d.mkdir(parents=True, exist_ok=True)

try:
    # STEP 2: Find efflux pump genes
    print("\n[STEP 2] Finding resistance genes with ABRicate...")
    finder = GeneFinder()
    
    # Process the E. coli K-12 genome
    genome_file = genomes_dir / "GCF_000005845.2.fasta"
    genes_df = finder.find_resistance_genes(fasta_file=genome_file)
    
    # Save results
    if not genes_df.empty:
        genes_df.to_csv(results_dir / "gene_finder.csv", index=False)
        print(f"✓ Found {len(genes_df)} resistance genes")
        print("\nGenes detected:")
        print(genes_df[['gene', 'product', 'coverage', 'identity']].to_string(index=False))
    else:
        print("⚠ No resistance genes found")
    
    # STEP 3: Extract and translate sequences
    print("\n[STEP 3] Extracting & translating sequences...")
    if not genes_df.empty:
        extractor = SequenceExtractor()
        for idx, row in genes_df.iterrows():
            try:
                protein_seq = extractor.extract_and_translate(
                    fasta_file=genome_file,
                    gene_name=row['gene'],
                    start=row['start'],
                    end=row['end'],
                    strand=row['strand']
                )
                
                # Write protein FASTA
                output_file = proteins_dir / f"{row['gene']}_{row['accession']}.faa"
                with open(output_file, 'w') as f:
                    f.write(f">{row['gene']}|{row['accession']}|{row['product']}\n")
                    f.write(str(protein_seq.seq) + "\n")
                    
            except Exception as e:
                print(f"  ⚠ Skipped {row['gene']}: {e}")
                continue
        
        protein_count = len(list(proteins_dir.glob("*.faa")))
        print(f"✓ Extracted {protein_count} protein sequences")
    else:
        print("⚠ No genes to extract")
    
    # STEP 4: Call variants with ML
    print("\n[STEP 4] Calling variants with ML predictions...")
    caller = VariantCaller(
        ml_model_path="/app/models/ciprofloxacin_predictor.pkl",
        enable_ml=True,
        antibiotic="Ciprofloxacin"
    )
    
    mutations_df = caller.call_variants(
        protein_dir=str(proteins_dir),
        reference_dir=str(refs_dir),
        output_csv=str(results_dir / "mutation_report.csv")
    )
    print(f"✓ Identified {len(mutations_df)} mutations")
    if not mutations_df.empty:
        print("\nMutation Report:")
        print(mutations_df[['gene', 'mutation', 'ml_prediction', 'ml_confidence']])
    
    # STEP 5: Visualize with PyMOL
    print("\n[STEP 5] Generating 3D visualizations...")
    visualizer = PyMOLVisualizer()
    vis_count = visualizer.visualize_mutations(
        mutations_df=mutations_df,
        output_dir=str(viz_dir)
    )
    print(f"✓ Generated {vis_count} visualizations")
    
    print("\n" + "="*70)
    print("  PIPELINE COMPLETE")
    print("="*70)
    print(f"\nOutputs:")
    print(f"  - Mutation report: {results_dir / 'mutation_report.csv'}")
    print(f"  - Visualizations:  {viz_dir}/")
    
except Exception as e:
    print(f"\n❌ Pipeline failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
