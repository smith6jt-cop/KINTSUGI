#!/usr/bin/env python3
"""Scan orange storage and build a dataset manifest from experiment.json files."""

import json, csv, sys
from pathlib import Path

ORANGE_ROOT = Path("/orange/maigan/HuBMAP_CODEX_LN-SP_Datasets")
OUTPUT_CSV = Path("/blue/maigan/smith6jt/dataset_manifest.csv")

# Possible names for the channel names file
CHANNEL_FILE_NAMES = [
    "channelnames.txt", "CHANNELNAMES.txt", "channel_names.txt",
    "channelnames.csv", "channels.txt", "markers.txt",
]

rows = []
errors = []
seen_data_roots = {}  # data_root.name -> dataset name (for dedup)

# Helper: try multiple field names, return first match or default
def get(exp, *keys, default="?"):
    for k in keys:
        v = exp.get(k)
        if v is not None:
            return v
    return default

all_dirs = sorted(d for d in ORANGE_ROOT.iterdir() if d.is_dir())
total = len(all_dirs)

for i, dataset_dir in enumerate(all_dirs, 1):
    # Look for experiment.json here or up to 3 levels deep
    exp_json = dataset_dir / "experiment.json"
    data_root = dataset_dir
    if not exp_json.exists():
        for depth in range(1, 4):
            matches = list(dataset_dir.glob("/".join(["*"] * depth + ["experiment.json"])))
            if matches:
                exp_json = matches[0]
                data_root = exp_json.parent
                break

    if not exp_json.exists():
        errors.append(f"SKIP  {dataset_dir.name}: no experiment.json")
        print(f"[{i}/{total}] SKIP {dataset_dir.name} (no experiment.json)")
        continue

    # Deduplicate: skip if we already processed a dataset with the same data_root name
    # (e.g., CX_19-002_lymph-node_R1/src_X and top-level src_X point to same data)
    root_key = data_root.name
    if root_key in seen_data_roots:
        prev = seen_data_roots[root_key]
        errors.append(f"SKIP  {dataset_dir.name}: duplicate of {prev} (same data_root {root_key})")
        print(f"[{i}/{total}] SKIP {dataset_dir.name} (duplicate of {prev})")
        continue
    seen_data_roots[root_key] = dataset_dir.name

    print(f"[{i}/{total}] {dataset_dir.name} ...", end=" ", flush=True)

    try:
        with open(exp_json) as f:
            exp = json.load(f)
    except json.JSONDecodeError as e:
        errors.append(f"SKIP  {dataset_dir.name}: bad JSON ({e.msg}, line {e.lineno})")
        print(f"bad JSON (line {e.lineno}: {e.msg})")
        continue

    # Find channel names file (check both top-level and data_root)
    ch_file = None
    for search_dir in dict.fromkeys([dataset_dir, data_root]):
        for name in CHANNEL_FILE_NAMES:
            candidate = search_dir / name
            if candidate.exists():
                ch_file = candidate.name
                break
        if ch_file:
            break

    # Count cycle directories: only dirs starting with "cyc"
    cycle_dirs = sorted([
        d.name for d in data_root.iterdir()
        if d.is_dir() and d.name.startswith("cyc")
    ])
    n_cycle_dirs = len(cycle_dirs)

    # Read parameters (try CODEX then KINTSUGI field names)
    n_cycles  = get(exp, "numCycles", "n_cycles", default=n_cycle_dirs)
    n_ch      = get(exp, "numChannels", "channels_per_cycle", default=4)
    tile_rows = get(exp, "regionHeight", "tile_rows", default=0)
    tile_cols = get(exp, "regionWidth", "tile_cols", default=0)
    n_zplanes = get(exp, "numZPlanes", "n_zplanes", default=0)

    # Flag cycle count mismatches
    if n_cycle_dirs != int(n_cycles) and n_cycle_dirs > 0:
        errors.append(f"WARN  {dataset_dir.name}: n_cycles={n_cycles} but {n_cycle_dirs} cyc dirs")

    # Estimate raw size: sample one tif, multiply by computed file count
    n_files = int(n_cycles) * int(n_ch) * int(tile_rows) * int(tile_cols) * int(n_zplanes)
    raw_gb = -1
    if n_files > 0 and cycle_dirs:
        try:
            sample = next((data_root / cycle_dirs[0]).rglob("*.tif"), None)
            if sample:
                raw_gb = (sample.stat().st_size * n_files) / 1e9
        except Exception:
            pass

    size_str = f"{raw_gb:.1f}GB" if raw_gb >= 0 else "?"
    print(f"cycles={n_cycles} ch={n_ch} grid={tile_rows}x{tile_cols} z={n_zplanes} est={size_str}")

    rows.append({
        "dataset":        dataset_dir.name,
        "orange_path":    str(data_root),
        "n_cycles":       n_cycles,
        "n_cycle_dirs":   n_cycle_dirs,
        "channels":       n_ch,
        "tile_rows":      tile_rows,
        "tile_cols":      tile_cols,
        "n_zplanes":      n_zplanes,
        "xy_pixel_nm":    get(exp, "xyResolution", "xy_pixel_size"),
        "z_step_nm":      get(exp, "zPitch", "z_step_size"),
        "NA":             get(exp, "aperture", "numerical_aperture"),
        "RI":             get(exp, "tissue_refractive_index"),
        "ch_names_file":  ch_file or "MISSING",
        "raw_size_GB":    f"{raw_gb:.1f}" if raw_gb >= 0 else "?",
        "first_cycle":    cycle_dirs[0] if cycle_dirs else "?",
    })

# Write CSV
fieldnames = list(rows[0].keys()) if rows else []
with open(OUTPUT_CSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(rows)

# Summary
total_gb = sum(float(r["raw_size_GB"]) for r in rows if r["raw_size_GB"] != "?")
print(f"\nManifest written: {OUTPUT_CSV}")
print(f"  {len(rows)} unique datasets")
print(f"  Total estimated raw size: {total_gb:.0f} GB ({total_gb/1000:.1f} TB)")
if errors:
    print(f"  {len(errors)} notes:")
    for e in errors:
        print(f"    {e}")
