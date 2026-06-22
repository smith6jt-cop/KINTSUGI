# Publishing this skill to the Skills Registry

This is a drafted `/retrospective` skill, preserved in the KINTSUGI repo because
the `Skills_Registry` is a **separate** repository
(`github.com/smith6jt-cop/Skills_Registry`). It already passed the registry's
validator locally (`OK mcp-fastmcp-structured-output`).

To publish it:

```bash
# from the KINTSUGI repo root
git submodule update --init Skills_Registry
cp -r docs/skills/mcp-fastmcp-structured-output Skills_Registry/plugins/kintsugi/

cd Skills_Registry
git checkout main && git pull            # leave the detached submodule HEAD
python scripts/generate_marketplace.py   # adds it to marketplace.json
python scripts/validate_plugins.py       # expect: OK mcp-fastmcp-structured-output
git add plugins/kintsugi/mcp-fastmcp-structured-output marketplace.json
git commit -m "feat(kintsugi): add mcp-fastmcp-structured-output skill"
git push
```

Optionally, bump the submodule pointer in KINTSUGI afterward:

```bash
cd ..                # back to KINTSUGI root
git add Skills_Registry
git commit -m "chore: bump Skills_Registry (add mcp-fastmcp-structured-output)"
```
