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

# Type checking
mypy src/ tests/
```

**Style guidelines:**
- Line length: 100 characters
- Formatter: Black
- Linter: Ruff (with isort, pyupgrade, flake8-bugbear)
- Python 3.10+ required

## Commit Message Format

This project uses [Conventional Commits](https://www.conventionalcommits.org/) for automatic semantic versioning.

```
<type>(<scope>): <description>

[optional body]

[optional footer(s)]
```

**Types:**
| Type | Description | Version Bump |
|------|-------------|-------------|
| `feat` | New feature | MINOR |
| `fix` | Bug fix | PATCH |
| `docs` | Documentation changes | - |
| `style` | Code style changes | - |
| `refactor` | Code refactoring | - |
| `perf` | Performance improvements | - |
| `test` | Adding/updating tests | - |
| `build` | Build system changes | - |
| `ci` | CI/CD changes | - |
| `chore` | Maintenance tasks | - |

**Breaking changes:** Add `!` after type or include `BREAKING CHANGE:` in footer (bumps MAJOR version).

**Examples:**
```bash
feat(mcp): add new image processing tool
fix(segmentation): resolve edge detection issue
docs: update installation instructions
refactor(signal)!: change API for background subtraction
```

## Pre-commit Hooks

Install hooks for automatic commit validation and code formatting:

```bash
pre-commit install
pre-commit install --hook-type commit-msg
```

This validates:
- Commit message format (Conventional Commits)
- Code formatting (Black)
- Linting (Ruff)

## Pull Request Process

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/my-feature`)
3. Make your changes
4. Run tests and linting
5. Commit using Conventional Commits format
6. Push to your fork (`git push origin feature/my-feature`)
7. Create a Pull Request

## Release Process

Releases are automated via GitHub Actions using Conventional Commits:

1. **Automatic releases**: Push to `main` triggers version analysis and release creation
2. **Manual releases**: Use workflow dispatch in GitHub Actions
3. **Local release script**:

```bash
# Preview release (dry run)
python scripts/release.py --dry-run --auto

# Manual bump
python scripts/release.py --bump minor
```

The `CHANGELOG.md` is automatically updated during releases.

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
