# MutationScan Examples

This directory contains example datasets, configuration files, and workflow demonstrations for the MutationScan antimicrobial resistance analysis pipeline.

## 📁 Contents

### Sample Data
- **`demo_accessions.txt`** - NCBI accession numbers for demonstration
- **`sample_genomes/`** - Small test genome files (when available)
- **`reference_sequences/`** - Wild-type reference sequences for alignment

### Configuration Examples
- **`basic_config.json`** - Basic pipeline configuration
- **`advanced_config.json`** - Advanced options and customization
- **`batch_processing.json`** - Configuration for large-scale analysis

### Workflow Examples
- **`quick_start_workflow.md`** - 15-minute getting started guide
- **`full_pipeline_example.md`** - Complete end-to-end analysis
- **`custom_analysis_example.md`** - Specialized analysis workflows

### Scripts
- **`run_demo.py`** - Automated demonstration script
- **`batch_analysis.py`** - Example batch processing script
- **`visualization_examples.py`** - Custom visualization scripts

## 🚀 Quick Start

### Basic Demo Run

```bash
# Navigate to the examples directory
cd examples

# Run a basic demonstration
python ../subscan/tools/run_harvester.py \
  --accessions-file demo_accessions.txt \
  --output-dir demo_output

# Continue with annotation
python ../subscan/tools/run_annotator.py \
  --manifest demo_output/genome_manifest.json \
  --output-dir demo_annotation
```

### Full Pipeline Example

```bash
# Run the complete pipeline demonstration
python run_demo.py --full-pipeline --output-dir complete_demo
```

## 📊 Expected Outputs

### Demo Dataset Statistics
- **Genomes**: 3-5 bacterial genomes
- **Expected Runtime**: 5-10 minutes
- **Output Size**: ~50MB
- **Mutations Expected**: 15-25 variants
- **Gene Families**: 3-5 resistance genes

### Example Results
- Interactive HTML report with mutation visualizations
- Gene co-occurrence analysis
- Multiple sequence alignments
- Publication-ready figures

## 🔧 Customization

### Using Your Own Data

1. **Replace accession file**:
   ```bash
   # Create your own accession list
   echo "your_accession_1" > my_accessions.txt
   echo "your_accession_2" >> my_accessions.txt
   ```

2. **Modify configuration**:
   ```json
   {
     "threads": 8,
     "quality_threshold": 0.95,
     "gene_families": ["gyrA", "parC", "rpoB"]
   }
   ```

3. **Run with custom settings**:
   ```bash
   python run_demo.py --config my_config.json --accessions my_accessions.txt
   ```

## 📚 Learning Resources

### Step-by-Step Tutorials
1. **[Basic Analysis](basic_analysis_tutorial.md)** - Simple resistance gene analysis
2. **[Advanced Features](advanced_features_tutorial.md)** - Custom workflows and optimization
3. **[Visualization Guide](visualization_tutorial.md)** - Creating publication-ready figures

### Best Practices
- Always validate input data quality
- Use appropriate reference sequences
- Monitor computational resources for large datasets
- Review results with domain expertise

## 🧪 Testing Examples

# Example Data & Demo Workflow Testing

To test the MutationScan pipeline with example data:

1. Create `examples/demo_accessions.txt` with a few NCBI accessions (one per line):
   - CP139196.1
   - CP178438.1
   - NZ_JAVRGK010000002.1
   - NZ_JBQGWU010000032.1
   - NZ_JBQGWV010000050.1

2. Create `examples/gene_list.txt` with a few gene names (one per line):
   - acrA
   - acrB
   - tolC

3. Run the pipeline:
   ```bash
   python subscan/tools/run_pipeline.py --accessions examples/demo_accessions.txt --gene-list examples/gene_list.txt --email you@example.com --output-dir data_output/demo_workflow --sepi-species "Escherichia coli" --threads 2 --verbose
   ```

4. Check for errors (logic, unicode, syntax, path). If any domino tool is missing, follow the README instructions to install.

5. After testing, delete the demo files to keep the repository clean:
   - `examples/demo_accessions.txt`
   - `examples/gene_list.txt`

You can repeat this process any time to validate the pipeline with sample data.

## 📞 Support

If you encounter issues with any examples:

1. Check the [Troubleshooting Guide](../docs/troubleshooting.md)
2. Review [GitHub Issues](https://github.com/vihaankulkarni29/MutationScan/issues)
3. Create a new issue with example details

---

**Happy analyzing!** 🧬