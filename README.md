# 🧬 MutationScan

**Analyze antimicrobial resistance mutations in bacterial genomes with just one command.**

[![Docker Hub](https://img.shields.io/badge/Docker%20Hub-mutationscan-2496ED?logo=docker)](https://hub.docker.com/r/vihaankulkarni29/mutationscan)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Version](https://img.shields.io/badge/version-1.0.0-green.svg)](https://github.com/vihaankulkarni29/MutationScan/releases)
[![Pipeline](https://img.shields.io/badge/pipeline-7%20stages-brightgreen.svg)](https://github.com/vihaankulkarni29/MutationScan)

MutationScan is a complete bioinformatics pipeline that automatically analyzes bacterial genomes for antimicrobial resistance mutations. Simply provide NCBI accession numbers and gene names—MutationScan handles the rest, generating beautiful interactive reports with mutation analysis, co-occurrence patterns, and publication-ready visualizations.

---

## 🚀 Quick Start

Get started in 3 simple steps:

### **Step 1: Install Docker Desktop**

Download and install Docker Desktop for your platform:

- **Windows**: [Download Docker Desktop](https://docs.docker.com/desktop/install/windows-install/)
- **macOS**: [Download Docker Desktop](https://docs.docker.com/desktop/install/mac-install/)
- **Linux**: [Install Docker Engine](https://docs.docker.com/engine/install/)

After installation, start Docker Desktop and ensure it's running.

---

### **Step 2: Prepare Your Data**

Create two folders on your computer and prepare your input files:

```bash
# Create folders
mkdir data_input data_output

# Navigate to data_input folder
cd data_input
```

Create two text files inside `data_input`:

**`accessions.txt`** - List your NCBI accession numbers (one per line):
```
NZ_CP107554
NZ_CP107555
NZ_CP107556
```

**`genes.txt`** - List target genes to analyze (one per line):
```
gyrA
parC
rpoB
```

---

### **Step 3: Run MutationScan**

Open your terminal (Command Prompt on Windows, Terminal on macOS/Linux) and run:

#### **Windows (PowerShell/CMD)**
```powershell
docker run -v C:\Users\YourName\data_input:/data/input -v C:\Users\YourName\data_output:/data/output vihaankulkarni29/mutationscan:latest
```

#### **macOS/Linux**
```bash
docker run -v ~/data_input:/data/input -v ~/data_output:/data/output vihaankulkarni29/mutationscan:latest
```

**Replace** `C:\Users\YourName\` or `~/` with the actual path to your folders.

---

### **� That's It!**

MutationScan will:
1. ✅ Download genomes from NCBI
2. ✅ Identify antimicrobial resistance genes
3. ✅ Extract and align sequences
4. ✅ Detect mutations (SNPs, insertions, deletions)
5. ✅ Analyze co-occurrence patterns
6. ✅ Generate interactive HTML reports

**Results will appear in your `data_output` folder**, including:
- 📊 Interactive HTML dashboard
- � Mutation analysis reports
- � CSV data files for further analysis
- 🧬 Sequence alignment files

---

## ✨ What MutationScan Does

MutationScan is a **7-stage automated pipeline** that transforms raw NCBI accession numbers into comprehensive AMR analysis:

| Stage | Name | Function |
|-------|------|----------|
| 🧬 **1** | Genome Harvester | Downloads genomes from NCBI with quality scoring |
| 🔍 **2** | Gene Annotator | Identifies AMR genes using CARD database |
| 🧪 **3** | Sequence Extractor | Extracts target gene sequences |
| 📏 **4** | Wild-type Aligner | Aligns sequences to reference genes |
| 🔬 **5** | Mutation Analyzer | Detects SNPs, insertions, and deletions |
| 📊 **6** | Co-occurrence Analyzer | Discovers resistance gene patterns |
| 📄 **7** | Report Generator | Creates interactive visualizations |

---

## 🎯 Key Features

### **For Researchers**
- 🔬 **Comprehensive Mutation Detection** - SNPs, insertions, deletions, and complex variants
- 📊 **Statistical Analysis** - Mutation frequencies, distributions, and significance testing
- 🧬 **Gene Co-occurrence Patterns** - Discover relationships between resistance genes
- 📈 **Publication-Ready Outputs** - High-quality figures and interactive dashboards

### **For Everyone**
- ⚡ **No Installation Required** - Everything runs in Docker
- 🖱️ **One Command** - Complete analysis from start to finish
- 🌐 **Works Everywhere** - Windows, macOS, and Linux
- 🎨 **Beautiful Reports** - Interactive HTML dashboards with Plotly.js visualizations

### **Technical Excellence**
- ✅ **Automated Quality Control** - Genome quality scoring and validation
- 🔄 **Reproducible Results** - Consistent analysis every time
- 🚀 **Fast Processing** - Multi-threaded analysis for large datasets
- 📚 **Comprehensive Data** - NCBI genomes + CARD AMR database

---

## 📊 Example Output

After running MutationScan, you'll receive:

### **Interactive Dashboard**
- 📈 Mutation distribution charts (interactive with Plotly.js)
- 🔗 Gene co-occurrence network diagrams
- 📋 Sortable/filterable data tables
- 🧬 Multiple sequence alignment viewer

### **Analysis Files**
- `mutation_analysis_report.html` - Main interactive dashboard
- `mutation_summary.csv` - All detected mutations
- `cooccurrence_matrix.csv` - Gene relationship data
- `genome_quality_scores.csv` - Quality metrics for each genome

### **Sample Results**
```
✓ Total Genomes Analyzed: 150
✓ Mutations Detected: 5,284
✓ Co-occurrence Pairs: 171,379
✓ Gene Families Analyzed: 12
```

---

## 💡 Tips & Best Practices

### **Input Preparation**
- ✅ Use valid NCBI accession numbers (RefSeq or GenBank format)
- ✅ One accession per line in `accessions.txt`
- ✅ Use standard gene names (e.g., gyrA, parC, rpoB)
- ✅ Check [CARD database](https://card.mcmaster.ca/) for supported AMR genes

### **Performance Optimization**
- 🚀 Analyze 10-50 genomes for fastest results (5-15 minutes)
- 🚀 100+ genomes may take 30-60 minutes depending on genome size
- 🚀 Use specific gene lists instead of "all genes" for faster analysis
- 🚀 Ensure stable internet connection for NCBI downloads

### **Troubleshooting**
- ❓ **Container exits immediately?** Check Docker Desktop is running
- ❓ **No results?** Verify input files are in `data_input` folder
- ❓ **Permission errors?** Ensure folders have read/write permissions
- ❓ **Network errors?** Check internet connection and NCBI API availability

---

## � Additional Resources

### **Documentation**
- 📖 **[Docker Usage Guide](DOCKER.md)** - Advanced Docker usage and troubleshooting
- 📖 **[Contributing Guidelines](CONTRIBUTING.md)** - For developers and contributors
- 📖 **[Changelog](CHANGELOG.md)** - Version history and release notes
- 📖 **[Security Policy](SECURITY.md)** - Security considerations

### **Example Datasets**
- 📂 **[Sample Data](examples/)** - Ready-to-use example files
- 📂 **[Tutorials](docs/)** - Step-by-step tutorials and guides

### **External Resources**
- 🌐 [NCBI Database](https://www.ncbi.nlm.nih.gov/) - Genome accession numbers
- 🌐 [CARD Database](https://card.mcmaster.ca/) - AMR gene information
- 🌐 [Docker Documentation](https://docs.docker.com/) - Docker help and tutorials

---

## 🆘 Getting Help

### **Need Support?**

1. **Check Documentation**: Review [DOCKER.md](DOCKER.md) for detailed troubleshooting
2. **Report Issues**: Found a bug? [Create an issue](https://github.com/vihaankulkarni29/MutationScan/issues/new/choose)
3. **Request Features**: Have an idea? [Submit a feature request](https://github.com/vihaankulkarni29/MutationScan/issues/new/choose)
4. **Ask Questions**: Join [GitHub Discussions](https://github.com/vihaankulkarni29/MutationScan/discussions)

### **Common Questions**

**Q: Do I need to know Python or bioinformatics?**  
A: No! MutationScan is designed for anyone. Just follow the 3-step Quick Start.

**Q: Can I analyze my own genome sequences?**  
A: Currently, MutationScan works with NCBI accession numbers. Local file support is planned.

**Q: How much does it cost?**  
A: MutationScan is completely free and open-source (MIT License).

**Q: Can I use this for publication?**  
A: Yes! See the Citation section below.

---

## 👨‍💻 For Developers

Want to contribute or customize MutationScan?

### **Development Setup**
See [CONTRIBUTING.md](CONTRIBUTING.md) for:
- Local development environment setup
- Code contribution guidelines
- Testing and quality standards
- Pull request process

### **Build from Source**
```bash
git clone https://github.com/vihaankulkarni29/MutationScan.git
cd MutationScan
docker build -t mutationscan:dev .
```

---

## 📜 License

MutationScan is open-source software licensed under the **MIT License**.

You are free to:
- ✅ Use commercially
- ✅ Modify and distribute
- ✅ Use privately
- ✅ Use in research and publications

See [LICENSE](LICENSE) file for full details.


## � Citation

If you use MutationScan in your research, please cite:

```bibtex
@software{kulkarni2025mutationscan,
  title={MutationScan: A Comprehensive Bioinformatics Pipeline for Antimicrobial Resistance Analysis},
  author={Kulkarni, Vihaan},
  year={2025},
  url={https://github.com/vihaankulkarni29/MutationScan},
  version={1.0.0},
  doi={10.5281/zenodo.XXXXXXX}
}
```

---

## 🙏 Acknowledgments

MutationScan stands on the shoulders of giants:

- **[NCBI](https://www.ncbi.nlm.nih.gov/)** - Comprehensive genomic databases
- **[CARD Database](https://card.mcmaster.ca/)** - Curated antimicrobial resistance gene data
- **[EMBOSS Suite](http://emboss.sourceforge.net/)** - Sequence alignment and analysis tools
- **[Plotly.js](https://plotly.com/javascript/)** - Interactive data visualizations
- **[Docker](https://www.docker.com/)** - Containerization platform

Thank you to all contributors and users making AMR research more accessible!

---

## 📧 Contact & Support

**Author**: Vihaan Kulkarni  
**Email**: [vihaankulkarni29@gmail.com](mailto:vihaankulkarni29@gmail.com)  
**GitHub**: [@vihaankulkarni29](https://github.com/vihaankulkarni29)

### Quick Links
- 🐛 [Report a Bug](https://github.com/vihaankulkarni29/MutationScan/issues/new?template=bug_report.yml)
- 💡 [Request a Feature](https://github.com/vihaankulkarni29/MutationScan/issues/new?template=feature_request.yml)
- 📖 [Improve Documentation](https://github.com/vihaankulkarni29/MutationScan/issues/new?template=documentation.yml)
- 💬 [Join Discussions](https://github.com/vihaankulkarni29/MutationScan/discussions)

---

<div align="center">

### 🧬 Advancing Antimicrobial Resistance Research 🧬

**MutationScan v1.0.0** | [Docker Hub](https://hub.docker.com/r/vihaankulkarni29/mutationscan) | [GitHub](https://github.com/vihaankulkarni29/MutationScan) | [Documentation](DOCKER.md)

*Made with ❤️ for the AMR research community*

</div>