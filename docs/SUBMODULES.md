# Git Submodules Guide

This document provides instructions for working with Git submodules in the KINTSUGI project.

## Current Submodules

| Submodule | Path | Purpose |
|-----------|------|---------|
| Skills_Registry | `Skills_Registry/` | Shared skills and learnings repository for Claude Code |

## Quick Start

### Cloning the Repository with Submodules

When cloning KINTSUGI for the first time, use one of these methods:

```bash
# Option 1: Clone with submodules in one command (recommended)
git clone --recurse-submodules https://github.com/smith6jt-cop/KINTSUGI.git

# Option 2: Clone then initialize submodules
git clone https://github.com/smith6jt-cop/KINTSUGI.git
cd KINTSUGI
git submodule update --init --recursive
```

### Pulling Updates with Submodules

Always include submodules when pulling:

```bash
git pull --recurse-submodules
```

## Common Operations

### Check Submodule Status

```bash
# View current submodule state
git submodule status

# For nested submodules (if any)
git submodule status --recursive
```

### Update Submodules to Tracked Commits

If teammates have updated submodule pointers:

```bash
# Sync submodules to commits tracked by parent repo
git submodule update
```

### Update Submodules to Latest Remote

To get the latest changes from submodule remotes:

```bash
# Update all submodules to their latest remote commits
git submodule update --remote

# Update a specific submodule
git submodule update --remote Skills_Registry
```

### Making Changes in a Submodule

1. Navigate to the submodule directory:
   ```bash
   cd Skills_Registry
   ```

2. Checkout the appropriate branch:
   ```bash
   git checkout main
   ```

3. Make your changes and commit:
   ```bash
   git add .
   git commit -m "feat: your changes"
   git push
   ```

4. Return to parent and commit the updated pointer:
   ```bash
   cd ..
   git add Skills_Registry
   git commit -m "chore: update Skills_Registry submodule"
   git push
   ```

## Recommended Git Aliases

Add these aliases to simplify submodule operations:

```bash
# Clone repositories with all submodules
git config --global alias.clone-all 'clone --recurse-submodules'

# Pull with automatic submodule updates
git config --global alias.pull-all 'pull --recurse-submodules'

# Update all submodules to tracked commits
git config --global alias.subup 'submodule update --init --recursive'

# Update all submodules to latest remote
git config --global alias.sublatest 'submodule update --remote --merge'
```

Usage:
```bash
git clone-all <repository-url>
git pull-all
git subup
git sublatest
```

## Troubleshooting

### Detached HEAD in Submodule

If `git status` shows `(modified content)` for a submodule:

```bash
git submodule update
```

### Submodule Directory is Empty

```bash
git submodule update --init --recursive
```

### Submodule Changes Not Reflected

After pulling, always run:

```bash
git submodule update
```

### Reset Submodule to Tracked Commit

If a submodule has unwanted local changes:

```bash
cd Skills_Registry
git checkout .
git clean -fd
cd ..
git submodule update
```

## Adding a New Submodule

```bash
# Basic addition
git submodule add <repository-url> <path>

# With specific branch tracking
git submodule add -b main <repository-url> <path>

# Commit the addition
git add .gitmodules <path>
git commit -m "chore: add <name> submodule"
```

## Removing a Submodule

```bash
# Remove submodule entry from .gitmodules
git config -f .gitmodules --remove-section submodule.<path>

# Remove submodule entry from .git/config
git config --remove-section submodule.<path>

# Remove submodule from index
git rm --cached <path>

# Remove submodule directory
rm -rf <path>
rm -rf .git/modules/<path>

# Commit the removal
git commit -m "chore: remove <name> submodule"
```

## Best Practices

1. **Always use `--recurse-submodules`** when cloning or pulling to ensure submodules stay in sync.

2. **Don't commit build artifacts** in submodules - they're system-specific and should be built locally.

3. **Update after teammate changes** - Run `git submodule update` after pulling when others have modified submodule pointers.

4. **Commit submodule pointer changes** - After updating a submodule, commit the updated reference in the parent repository.

5. **Use branch tracking** for actively developed submodules:
   ```bash
   git config -f .gitmodules submodule.Skills_Registry.branch main
   ```

6. **Keep submodules shallow** if history isn't needed:
   ```bash
   git submodule update --init --depth 1
   ```

## References

- [Git Submodules Documentation](https://git-scm.com/book/en/v2/Git-Tools-Submodules)
- [Git Submodules Basic Explanation](https://gist.github.com/gitaarik/8735255)
