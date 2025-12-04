# Contributing to KINTSUGI

Thank you for your interest in contributing to KINTSUGI!

## Development Setup

```bash
git clone https://github.com/smith6jt-cop/KINTSUGI.git
cd KINTSUGI
conda env create -f envs/env-linux.yml
conda activate KINTSUGI
pip install -e ".[dev]"
```

## Running Tests

```bash
# Run all tests
pytest tests/ -v

# Run with coverage
pytest tests/ --cov=src/kintsugi --cov-report=html
```

## Code Quality

```bash
# Lint code
ruff check src/ tests/

# Format code
black src/ tests/
```

## Pull Request Process

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/my-feature`)
3. Make your changes
4. Run tests and linting
5. Commit your changes (`git commit -am 'Add new feature'`)
6. Push to your fork (`git push origin feature/my-feature`)
7. Create a Pull Request

## Documentation

Documentation is built with Sphinx. To build locally:

```bash
cd docs
pip install -r requirements.txt

# On Linux/macOS
make html

# On Windows
.\make.bat html
```

The built documentation will be in `docs/_build/html/`.

## Reporting Issues

Please report issues at: https://github.com/smith6jt-cop/KINTSUGI/issues

Include:
- Operating system and version
- Python version
- Full error traceback
- Steps to reproduce
- Output of `kintsugi check`
