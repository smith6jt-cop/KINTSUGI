# Documentation Maintenance Guide

This guide explains how to build, maintain, and deploy the KINTSUGI documentation.

## Documentation Structure

```
docs/
├── conf.py              # Sphinx configuration
├── index.rst            # Main documentation page
├── Makefile             # Unix build script
├── make.bat             # Windows build script
├── requirements.txt     # Documentation dependencies
├── _static/             # Static files (CSS, images)
├── _templates/          # Custom templates
├── _build/              # Build output (gitignored)
├── installation.md      # Installation guide
├── quickstart.md        # Quick start guide
├── workflows.md         # Processing workflows
├── cli.md               # CLI reference
├── api.md               # API reference
├── contributing.md      # Contribution guide
├── TROUBLESHOOTING.md   # Troubleshooting guide
├── DEPENDENCY_ANALYSIS.md # Dependency analysis
├── DOCS_MAINTENANCE.md  # This file
├── CD3e.gif             # Demo animation
└── CD31_curtain.gif     # Demo animation
```

## Building Documentation Locally

### Prerequisites

Install documentation dependencies:

```bash
cd docs
pip install -r requirements.txt
```

### On Windows

**Option 1: Using make.bat**

```batch
cd docs
.\make.bat html
```

**Option 2: Direct sphinx-build command**

If `make` is not available (common on Windows), use sphinx-build directly:

```batch
cd docs
sphinx-build -M html . _build
```

**Option 3: Using Python module**

```batch
cd docs
python -m sphinx -M html . _build
```

### On Linux/macOS

```bash
cd docs
make html
```

### Viewing the Built Documentation

After building, open `docs/_build/html/index.html` in your web browser.

## Live Preview During Development

For live-reloading during documentation development:

```bash
cd docs
sphinx-autobuild . _build/html
```

Then open http://127.0.0.1:8000 in your browser. Changes will auto-reload.

## Deploying to GitHub Pages

### Automatic Deployment

Documentation is automatically deployed to GitHub Pages when changes are pushed to the main branch. The deployment is handled by the `.github/workflows/docs.yml` workflow.

### Manual Deployment

If needed, you can manually trigger deployment:

1. Go to the repository on GitHub
2. Navigate to Actions > Deploy Documentation
3. Click "Run workflow"

### Accessing the Live Site

Once deployed, the documentation is available at:

```
https://smith6jt-cop.github.io/KINTSUGI/
```

### Enabling GitHub Pages (First Time Setup)

1. Go to repository Settings > Pages
2. Under "Build and deployment", select:
   - Source: "GitHub Actions"
3. The workflow will handle the rest

## Adding New Documentation

### Adding a New Page

1. Create a new `.md` or `.rst` file in the `docs/` directory
2. Add the file to the appropriate `toctree` in `index.rst`

Example for `docs/new-page.md`:

```rst
.. toctree::
   :maxdepth: 2
   :caption: User Guide

   workflows
   new-page    # Add your new page here
```

### Adding API Documentation

For automatic API documentation from docstrings:

1. Ensure modules have proper docstrings (Google or NumPy style)
2. Use autodoc directives in `.rst` files:

```rst
.. automodule:: kintsugi.kreg
   :members:
   :undoc-members:
```

### Adding Images

1. Place images in `docs/_static/` or the `docs/` directory
2. Reference in Markdown:

```markdown
![Alt text](image.png)
```

Or in reStructuredText:

```rst
.. image:: image.png
   :alt: Alt text
```

## Writing Style Guidelines

1. Use clear, concise language
2. Include code examples where helpful
3. Keep paragraphs short
4. Use headers to organize content
5. Link to related documentation sections
6. Test all code examples before committing

## Common Tasks

### Clean Build

To force a complete rebuild:

```bash
# Unix
make clean && make html

# Windows
.\make.bat clean
.\make.bat html
```

### Check for Broken Links

```bash
make linkcheck
```

### Generate PDF (requires LaTeX)

```bash
make latexpdf
```

## Troubleshooting Documentation Builds

### "sphinx-build not found"

Install Sphinx:

```bash
pip install sphinx sphinx-rtd-theme myst-parser
```

### "make is not recognized" (Windows)

Use the direct sphinx-build command:

```batch
sphinx-build -M html . _build
```

Or install GNU Make via:
- [Chocolatey](https://chocolatey.org/): `choco install make`
- [Git for Windows](https://gitforwindows.org/) (includes make in Git Bash)

### Import errors during build

Ensure you're in the KINTSUGI conda environment:

```bash
conda activate KINTSUGI
cd docs
make html
```

### Theme not found

```bash
pip install sphinx-rtd-theme
```

### MyST parser errors

```bash
pip install myst-parser
```

## Version Updates

When releasing a new version:

1. Update `release` in `docs/conf.py`
2. Update `__version__` in `src/kintsugi/__init__.py`
3. Add release notes to documentation if applicable

## Contact

For documentation issues, open an issue at:
https://github.com/smith6jt-cop/KINTSUGI/issues
