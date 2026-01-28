# Documentation Directory

This directory contains comprehensive documentation for MutationScan.

## Contents

- **API_REFERENCE.md** - Complete API documentation for all modules
- **INSTALLATION.md** - Detailed installation instructions
- **QUICKSTART.md** - Quick start guide for new users
- **ARCHITECTURE.md** - System architecture and design patterns
- **TUTORIALS/** - Step-by-step tutorials for common workflows
- **EXAMPLES/** - Code examples and sample workflows

## Building Documentation

To build HTML documentation:

```bash
cd docs
pip install -e ".[docs]"
make html
```

Documentation will be available in `_build/html/index.html`

## Contributing to Documentation

- Update docstrings in source code (Google style)
- Add examples in module documentation
- Update tutorials for new features
- Keep API reference up to date
