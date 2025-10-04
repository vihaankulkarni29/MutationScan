# MutationScan Advanced Documentation

This file contains all advanced, technical, and developer-facing documentation for MutationScan. For user setup, pipeline usage, and troubleshooting, see the main README.md.

---

## ⚙️ Advanced Configuration

- All pipeline parameters can be set via CLI flags (see main README).
- For custom workflows, modify input files or domino tool arguments as needed.
- Environment variables for advanced users:
  - `MUTATIONSCAN_THREADS`: Override default thread count
  - `MUTATIONSCAN_OUTPUT_DIR`: Set default output directory

---

## 👨‍💻 Developer Guide

- Clone the repo and install with `pip install -e .[dev]`
- All domino tools are in `subscan/tools/` and `subscan/ncbi_genome_extractor/`
- To contribute:
  1. Fork the repo
  2. Create a feature branch
  3. Submit a pull request
- Code style: Black, Flake8, Mypy
- Tests: Pytest (see `subscan/tests/`)

---

## 🧬 API Reference

- Main entry point: `subscan.tools.run_pipeline:main`
- Each domino tool exposes a CLI and Python API (see source code for details)
- Manifest files are JSON, CSV, and HTML (see pipeline overview in main README)

---

## 🚀 Performance Tuning

- Use multi-threading for large datasets (`--threads` flag)
- Analyze 10-50 genomes for fastest results
- Ensure SSD storage and stable internet for optimal performance

---

## 🤝 Contribution & Security

- All contributions welcome via GitHub pull requests
- Security issues: Please report via [GitHub Issues](https://github.com/vihaankulkarni29/MutationScan/issues)
- MIT License: Free for academic, commercial, and personal use

---

## 📚 Additional Resources

- Main README: [../README.md](../README.md)
- Example Data: [../examples/](../examples/)
- Contact: vihaankulkarni29@gmail.com

---

*This file replaces all other technical markdowns. For any topic not covered, open a GitHub Issue or discussion.*