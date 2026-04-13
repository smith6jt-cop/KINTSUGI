#!/usr/bin/env python3
"""Stage HuBMAP CODEX datasets from PATH lab share to HiPerGator via Globus.

Source: PATH lab SMB share (path.ahc.ufl.edu) via Globus Connect Personal
Destination: HiPerGator /blue storage

The PATH share organizes datasets under release folders:
  Release1/  → Cases 1-4  (CX_19-001 through CX_19-004)
  Release2/  → Cases 5-9  (CX_20-005 through CX_20-009)
  Release3/  → Cases 10-14 (CX_21-010 through CX_21-015)

Each dataset has a src_* subdirectory containing cycle dirs + experiment.json.
This script maps manifest entries to the correct nested paths on the source.

Usage:
    module load globus
    globus login   # first time only

    python stage_datasets_globus.py list              # show all datasets & status
    python stage_datasets_globus.py transfer --all     # transfer all unstaged
    python stage_datasets_globus.py transfer --batch 5 # 5 smallest unstaged
    python stage_datasets_globus.py transfer --datasets CX_19-003_spleen_CC3-E
    python stage_datasets_globus.py transfer --all --dry-run
    python stage_datasets_globus.py status             # check active transfers
"""

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path, PurePosixPath

# ── Endpoints ──────────────────────────────────────────────────────────────────
SRC_ENDPOINT = "f1b69b9e-f07a-11ef-8c40-0e26ca329435"
# Use cifs mount, NOT GVFS — GVFS causes EOF errors / Path-not-allowed under load
SRC_BASE = "/mnt/ahc_share/SHARE/HuBMAP/"

DST_ENDPOINT = "5dbaf795-8a7e-4dca-91aa-6e10d610c2b3"
DST_BASE = "/blue/maigan/smith6jt/KINTSUGI_Projects/"

# Release folder names on the PATH share
RELEASE1 = "Codex_dataset_hubmap_Case1-4_Release1_D2019-2020"
RELEASE2 = "CODEX_dataset_hubmap_Case5-9_Release2"
RELEASE3 = "Codex_dataset_HuBMAP_cases10-14_release3"

# ── Source path mapping ────────────────────────────────────────────────────────
# Maps manifest dataset name → relative path from SRC_BASE to the src_* directory
# containing cycle dirs. Built from manual inspection of the PATH share structure.
SOURCE_PATHS = {
    # ── Release 1 (Cases 1-4) ──────────────────────────────────────────────
    "CX_19-001_SP_CC2-A28":
        f"{RELEASE1}/CX_19-001_SP_CC2-A28/data/src_CX_19-001-SP_CC2-A28",
    "CX_19-001_spleen_CC2-A":
        f"{RELEASE1}/CX_19-001_spleen_CC2-A/src_CX_19-001_spleen_CC2-A_D200211",
    "CX_19-001_spleen_CC3-C":
        f"{RELEASE1}/CX_19-001_spleen_CC3-C/src_CX_19-001_spleen_CC3-C",
    "CX_19-002_lymph-node_R1":
        f"{RELEASE1}/CX_19-002_lymph-node_R1/src_CX_19-002_lymph-node_R1",
    "CX_19-002_lymph-node_R2":
        f"{RELEASE1}/CX_19-002_lymph-node_R2/src_CX_19-002_lymph-node_R2",
    "CX_19-002_lymph-node_R3":
        f"{RELEASE1}/CX_19-002_lymph-node_R3/src_CX_19-002_lymph-node_R3",
    "CX_19-003_SP_CC2-E28":
        f"{RELEASE1}/CX_19-003_SP_CC2-E28/data/src_CX_19-003_SP_CC2-E28",
    "CX_19-003_lymph-node_R1":
        f"{RELEASE1}/CX_19-003_lymph-node_R1/src_CX_19-003_LN_1",
    "CX_19-003_lymph-node_R2":
        f"{RELEASE1}/CX_19-003_lymph-node_R2/src_CX_19-003_LN_R2",
    "CX_19-003_lymph-node_R3":
        f"{RELEASE1}/CX_19-003_lymph-node_R3/src_CX_19-003_LN_R3",
    "CX_19-003_spleen_CC1-A":
        f"{RELEASE1}/CX_19-003_spleen_CC1-A/src_CX_19-003_spleen_CC1-A",
    "CX_19-003_spleen_CC2-E":
        f"{RELEASE1}/CX_19-003_spleen_CC2-E/src_CX_19-003_spleen_CC2-E",
    "CX_19-003_spleen_CC3-E":
        f"{RELEASE1}/CX_19-003_spleen_CC3-E/src_CX_19-003_spleen_CC3-E_NEW",
    "CX_19-004_lymph-node_R1":
        f"{RELEASE1}/CX_19-004_lymph-node_R1/src_CX_19-004_LN_R1",
    "src_CX_19-002_spleen_CC2-A_D200210":
        f"{RELEASE1}/CX_19-002_spleen_CC2-A/src_CX_19-002_spleen_CC2-A_D200210",
    "src_CX_19-002_spleen_CC3-C":
        f"{RELEASE1}/CX_19-002_spleen_CC3-C/src_CX_19-002_spleen_CC3-C",
    "src_CX_19-004_LN_R2":
        f"{RELEASE1}/CX_19-004_lymph-node_R2/src_CX_19-004_LN_R2",
    "src_CX_19-004_LN_R3":
        f"{RELEASE1}/CX_19-004_lymph-node_R3/src_CX_19-004_LN_R3",
    "src_CX_19-004_spleen_CC1-1L":
        f"{RELEASE1}/CX_19-004_spleen_CC1-1L/src_CX_19-004_spleen_CC1-1L",
    "src_CX_19-004_spleen_CC2-2L":
        f"{RELEASE1}/CX_19-004_spleen_CC2-2L/src_CX_19-004_spleen_CC2-2L",
    "src_CX_19-004_spleen_CC3-2L":
        f"{RELEASE1}/CX_19-004_spleen_CC3-2L/src_CX_19-004_spleen_CC3-2L",
    # ── Release 2 (Cases 5-9) ─────────────────────────────────────────────
    "CX_20-008_LN_n10_reg001":
        f"{RELEASE2}/CX_20-008_LN/CX_20-008_LN_n10_reg001/src_CX_20-008_lymphnode_n10_reg001",
    "CX_20-008_LN_n9_reg002":
        f"{RELEASE2}/CX_20-008_LN/CX_20-008_LN_n9_reg002/src_CX_20-008_lymphnode_n9_reg002",
    "src_CX_20-007_SP_CC1-A":
        f"{RELEASE2}/CX_20-007_SP_CC1-A/src_CX_20-007_SP_CC1-A",
    "src_CX_20-007_SP_CC2-B":
        f"{RELEASE2}/CX_20-007_SP_CC2-B/src_CX_20-007_SP_CC2-B",
    "src_CX_20-008_SP_CC3-A":
        f"{RELEASE2}/CX_20-008_SP_CC3-A/src_CX_20-008_SP_CC3-A",
    # ── Release 3 (Cases 10-14) ───────────────────────────────────────────
    "src_CX_21-010_LN_n3":
        f"{RELEASE3}/CX_21-010_LN_n3/data/src_CX_21-010_LN_n3",
    "src_CX_21-010_SP-CC1-A":
        f"{RELEASE3}/CX_21-010_SP-CC1-A/src_CX_21-010_SP-CC1-A",
    "src_CX_21-012_LN_n3":
        f"{RELEASE3}/CX_21-012_LN_n3/src_CX_21-012_LN_n3",
    "src_CX_21-012_SP-CC3-A":
        f"{RELEASE3}/CX_21-012_SP-CC3-A/src_CX_21-012_SP-CC3-A",
    "src_CX_21-014_LN_2.5":
        f"{RELEASE3}/CX_21-014_LN_2.5/src_CX_21-014_LN_2.5",
    "src_CX_21-014_SP-CC2-A":
        f"{RELEASE3}/CX_21-014_SP-CC2-A/src_CX_21-014_SP-CC2-A",
    "src_CX_21-015_LN_n3_n1":
        f"{RELEASE3}/CX_21-015_LN_n3_n1/src_CX_21-015_LN_n3_n1",
    "src_CX_21-015_SP-CC3-B":
        f"{RELEASE3}/CX_21-015_SP-CC3-B/src_CX_21-015_SP-CC3-B",
    # ── Thymus datasets ───────────────────────────────────────────────────
    # Release 1
    "CX_19-001_thymus_CC1-A":
        f"{RELEASE1}/CX_19-001_thymus_CC1-A/src_CX_19-001_thymus_CC1-A",
    "CX_19-001_thymus_CC2-C":
        f"{RELEASE1}/CX_19-001_thymus_CC2-C/src_CX_19-001_thymus_CC2-C_NBF",
    "CX_19-002_thymus_CC1-A":
        f"{RELEASE1}/CX_19-002_thymus_CC1-A/src_CX_19-002_thymus_CC1-A",
    "CX_19-002_thymus_CC1-E":
        f"{RELEASE1}/CX_19-002_thymus_CC1-E/src_CX_19-002_thymus_CC1-E",
    "CX_19-002_thymus_CC2-C":
        f"{RELEASE1}/CX_19-002_thymus_CC2-C/src_CX_19-002_thymus_CC2-C",
    # Release 2
    "CX_20-005_TH_CC1-A":
        f"{RELEASE2}/CX_20-005_TH_CC1-A/src_CX_20-005_TH_CC1-A",
    "CX_20-005_TH_CC2-B":
        f"{RELEASE2}/CX_20-005_TH_CC2-B/src_CX_20-005_TH_CC2-B",
    "CX_20-006_TH_CC2-B":
        f"{RELEASE2}/CX_20-006_TH_CC2-B(NEW)/src_CX_20-006_TH_CC2-B",
    "CX_20-006_TH_CC2-C":
        f"{RELEASE2}/CX_20-006_TH_CC2-C/src_CX_20-006_TH_CC2-C",
    "CX_20-008_TH_CC1-A":
        f"{RELEASE2}/CX_20-008_TH_CC1-A/src_CX_20-008_TH_CC1-A",
    "CX_20-008_TH_CC2-B":
        f"{RELEASE2}/CX_20-008_TH_CC2-B/src_CX_20-008_TH_CC2-B",
    # Release 3
    "CX_21-011_TH_L-A":
        f"{RELEASE3}/CX_21-011_TH_L-A/src_CX_21-011_TH_L-A",
    "CX_21-011_L-B_TH_reg1-2":
        f"{RELEASE3}/CX_21-011 L-B_TH_reg1-2/src_CX_21-011 L-B_TH_reg1-2",
}

# ── Local paths ────────────────────────────────────────────────────────────────
MANIFEST = Path("/blue/maigan/smith6jt/dataset_manifest.csv")
PROJECTS_DIR = Path("/blue/maigan/smith6jt/KINTSUGI_Projects")
TASK_LOG = Path("/blue/maigan/smith6jt/globus_transfers.json")


def run_globus(*args, check=True):
    """Run a globus CLI command, return stdout."""
    cmd = ["globus"] + list(args)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if check and result.returncode != 0:
        stderr = result.stderr.strip()
        if "MissingLoginError" in stderr:
            print("ERROR: Not logged in. Run: module load globus && globus login")
            sys.exit(1)
        print(f"ERROR: globus {' '.join(args[:2])}...\n{stderr}", file=sys.stderr)
        sys.exit(1)
    return result.stdout.strip()


def load_manifest():
    """Load dataset manifest CSV into a dict keyed by dataset name."""
    if not MANIFEST.exists():
        print(f"WARNING: Manifest not found: {MANIFEST}")
        return {}
    datasets = {}
    with open(MANIFEST) as f:
        for row in csv.DictReader(f):
            datasets[row["dataset"]] = row
    return datasets


def get_staged_datasets():
    """Return set of dataset names that already have raw data."""
    staged = set()
    if not PROJECTS_DIR.exists():
        return staged
    for d in PROJECTS_DIR.iterdir():
        if not d.is_dir():
            continue
        marker = d / "data" / "raw" / ".staged"
        if marker.exists():
            staged.add(d.name)
    return staged


def load_task_log():
    """Load the Globus transfer task log."""
    if TASK_LOG.exists():
        with open(TASK_LOG) as f:
            return json.load(f)
    return {}


def save_task_log(log):
    """Save the Globus transfer task log."""
    with open(TASK_LOG, "w") as f:
        json.dump(log, f, indent=2)


def get_source_path(dataset_name):
    """Get the full Globus source path for a dataset."""
    rel = SOURCE_PATHS.get(dataset_name)
    if rel is None:
        return None
    return f"{SRC_BASE}{rel}/"


# ── Commands ───────────────────────────────────────────────────────────────────


def cmd_list(args):
    """List all manifest datasets and their staging status."""
    manifest = load_manifest()
    staged = get_staged_datasets()
    task_log = load_task_log()

    print(f"{'Dataset':<45} {'Size GB':<10} {'Staged':<10} {'Transfer':<15} {'Source Path'}")
    print("-" * 130)

    n_mapped = 0
    n_unmapped = 0
    for name, info in sorted(manifest.items()):
        size = info.get("raw_size_GB", "?")
        is_staged = "STAGED" if name in staged else "-"
        transfer = task_log.get(name, {}).get("status", "-")
        src_rel = SOURCE_PATHS.get(name)
        if src_rel:
            n_mapped += 1
            src_short = src_rel.split("/")[-1]  # just show src_* dir name
        else:
            n_unmapped += 1
            src_short = "NOT MAPPED"
        print(f"{name:<45} {size:<10} {is_staged:<10} {transfer:<15} {src_short}")

    print(f"\nTotal: {len(manifest)} datasets")
    print(f"  Mapped to source: {n_mapped}")
    print(f"  Not mapped: {n_unmapped}")
    print(f"  Already staged: {len(staged)}")
    unstaged = [n for n in manifest if n not in staged and n in SOURCE_PATHS]
    print(f"  Ready to transfer: {len(unstaged)}")


def cmd_transfer(args):
    """Submit Globus transfer tasks for datasets."""
    manifest = load_manifest()
    staged = get_staged_datasets()
    task_log = load_task_log()

    # Determine which datasets to transfer
    if args.datasets:
        to_transfer = args.datasets
    elif args.all or args.batch:
        # Get all unstaged, mapped datasets not already in task log
        unstaged = [
            name for name in manifest
            if name not in staged
            and name not in task_log
            and name in SOURCE_PATHS
        ]
        # Sort by size (smallest first) for batch mode
        unstaged.sort(
            key=lambda n: float(manifest[n].get("raw_size_GB", "9999"))
        )
        if args.batch:
            to_transfer = unstaged[: args.batch]
        else:
            to_transfer = unstaged
    else:
        print("ERROR: Specify --datasets, --all, or --batch N")
        sys.exit(1)

    if not to_transfer:
        print("Nothing to transfer — all datasets are staged or in progress.")
        return

    # Validate all datasets have source paths
    missing = [n for n in to_transfer if n not in SOURCE_PATHS]
    if missing:
        print(f"ERROR: No source path mapping for: {', '.join(missing)}")
        print("Add entries to SOURCE_PATHS dict in the script.")
        sys.exit(1)

    # Show plan
    total_gb = 0
    print(f"\n{'#':<4} {'Dataset':<45} {'Size GB':<10} {'Source (src_* dir)'}")
    print("-" * 100)
    for i, name in enumerate(to_transfer, 1):
        size = manifest.get(name, {}).get("raw_size_GB", "?")
        if size != "?":
            total_gb += float(size)
        src_rel = SOURCE_PATHS[name]
        src_short = PurePosixPath(src_rel).name
        print(f"{i:<4} {name:<45} {size:<10} {src_short}")
    print(f"\nTotal: {len(to_transfer)} datasets, ~{total_gb:.0f} GB")
    print(f"Destination: {DST_ENDPOINT}:{DST_BASE}{{dataset}}/data/raw/")

    if args.dry_run:
        print("\n[DRY RUN] No transfers submitted.")
        return

    # Confirm
    if not args.yes:
        resp = input("\nSubmit transfers? [y/N] ").strip().lower()
        if resp != "y":
            print("Aborted.")
            return

    # Submit one transfer task per dataset (allows independent tracking)
    submitted = 0
    for name in to_transfer:
        src_path = get_source_path(name)
        dst_path = f"{DST_BASE}{name}/data/raw/"

        print(f"\nTransferring: {name}")
        print(f"  src: {SRC_ENDPOINT}:{src_path}")
        print(f"  dst: {DST_ENDPOINT}:{dst_path}")

        try:
            output = run_globus(
                "transfer",
                f"{SRC_ENDPOINT}:{src_path}",
                f"{DST_ENDPOINT}:{dst_path}",
                "--recursive",
                "--label", f"KINTSUGI stage {name}",
                "--sync-level", "size",
                "--preserve-timestamp",
                "--notify", "failed,inactive",
            )
            # Parse task ID from output
            task_id = None
            for line in output.splitlines():
                if "Task ID:" in line:
                    task_id = line.split("Task ID:")[-1].strip()
            if not task_id:
                task_id = output.strip().split("\n")[-1].strip()

            print(f"  Task ID: {task_id}")
            submitted += 1

            task_log[name] = {
                "task_id": task_id,
                "status": "SUBMITTED",
                "src_path": src_path,
                "dst_path": dst_path,
                "size_gb": manifest.get(name, {}).get("raw_size_GB", "?"),
            }
            save_task_log(task_log)

        except SystemExit:
            print(f"  FAILED to submit transfer for {name}")
            task_log[name] = {"status": "FAILED_SUBMIT"}
            save_task_log(task_log)

    print(f"\nSubmitted {submitted}/{len(to_transfer)} transfer tasks.")
    print(f"Task log: {TASK_LOG}")
    print("Check status: python stage_datasets_globus.py status")


def cmd_status(args):
    """Check status of all Globus transfer tasks."""
    task_log = load_task_log()
    staged = get_staged_datasets()

    if not task_log:
        print("No transfer tasks recorded.")
        print(f"(Task log: {TASK_LOG})")
        return

    print(f"{'Dataset':<45} {'Task ID':<38} {'Status':<20} {'Size GB'}")
    print("-" * 115)

    active = 0
    completed = 0
    failed = 0

    for name, info in sorted(task_log.items()):
        task_id = info.get("task_id", "")
        size = info.get("size_gb", "?")

        if name in staged:
            status = "STAGED"
            completed += 1
        elif task_id and info.get("status") != "FAILED_SUBMIT":
            try:
                output = run_globus(
                    "task", "show", task_id,
                    "--format", "json",
                    check=False,
                )
                if output:
                    task_info = json.loads(output)
                    status = task_info.get("status", "UNKNOWN")
                    if status == "SUCCEEDED":
                        completed += 1
                        # Auto-create .staged marker
                        marker = PROJECTS_DIR / name / "data" / "raw" / ".staged"
                        if not marker.exists():
                            marker.parent.mkdir(parents=True, exist_ok=True)
                            marker.write_text(
                                f"staged via globus from {SRC_ENDPOINT}:{info.get('src_path', '')}\n"
                                f"task_id: {task_id}\n"
                            )
                    elif status == "ACTIVE":
                        active += 1
                        bt = task_info.get("bytes_transferred", 0)
                        bt_gb = bt / 1e9
                        status = f"ACTIVE ({bt_gb:.1f} GB)"
                    elif status in ("FAILED", "INACTIVE"):
                        failed += 1
                else:
                    status = info.get("status", "UNKNOWN")
            except (json.JSONDecodeError, Exception):
                status = info.get("status", "UNKNOWN")
        else:
            status = info.get("status", "UNKNOWN")
            if "FAIL" in status:
                failed += 1

        print(f"{name:<45} {task_id:<38} {status:<20} {size}")

    print(f"\nTotal: {len(task_log)} tasks")
    print(f"  Active: {active}  Completed: {completed}  Failed: {failed}")

    if failed > 0:
        print("\nTo retry failed: reset then re-transfer:")
        print("  python stage_datasets_globus.py reset --failed")
        print("  python stage_datasets_globus.py transfer --all")


def cmd_cancel(args):
    """Cancel active transfer tasks."""
    task_log = load_task_log()
    targets = args.datasets if args.datasets else list(task_log.keys())
    cancelled = 0

    for name in targets:
        info = task_log.get(name, {})
        task_id = info.get("task_id")
        if not task_id:
            continue
        print(f"Cancelling: {name} (task {task_id})")
        try:
            run_globus("task", "cancel", task_id, check=False)
            task_log[name]["status"] = "CANCELLED"
            cancelled += 1
        except Exception as e:
            print(f"  Failed: {e}")

    save_task_log(task_log)
    print(f"Cancelled {cancelled} tasks.")


def cmd_reset(args):
    """Clear task log entries for retrying."""
    task_log = load_task_log()

    if args.datasets:
        targets = args.datasets
    elif args.failed:
        targets = [
            name for name, info in task_log.items()
            if "FAIL" in info.get("status", "")
        ]
    elif args.all_entries:
        targets = list(task_log.keys())
    else:
        print("Specify --datasets, --failed, or --all")
        return

    for name in targets:
        task_log.pop(name, None)
        print(f"Cleared: {name}")

    save_task_log(task_log)


def main():
    parser = argparse.ArgumentParser(
        description="Stage HuBMAP datasets via Globus transfer"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # list
    sub.add_parser("list", help="List datasets and their status")

    # transfer
    p_xfer = sub.add_parser("transfer", help="Submit Globus transfer tasks")
    p_xfer.add_argument("--datasets", nargs="+", help="Specific dataset names")
    p_xfer.add_argument("--all", action="store_true", help="Transfer all unstaged")
    p_xfer.add_argument("--batch", type=int, metavar="N",
                        help="Transfer N smallest unstaged datasets")
    p_xfer.add_argument("--dry-run", action="store_true", help="Show plan only")
    p_xfer.add_argument("-y", "--yes", action="store_true", help="Skip confirmation")

    # status
    sub.add_parser("status", help="Check transfer task status")

    # cancel
    p_cancel = sub.add_parser("cancel", help="Cancel active transfers")
    p_cancel.add_argument("--datasets", nargs="+", help="Specific datasets (default: all)")

    # reset
    p_reset = sub.add_parser("reset", help="Clear task log entries for retry")
    p_reset.add_argument("--datasets", nargs="+", help="Specific datasets to clear")
    p_reset.add_argument("--failed", action="store_true", help="Clear all failed entries")
    p_reset.add_argument("--all", dest="all_entries", action="store_true",
                         help="Clear all entries")

    args = parser.parse_args()
    cmd = {
        "list": cmd_list,
        "transfer": cmd_transfer,
        "status": cmd_status,
        "cancel": cmd_cancel,
        "reset": cmd_reset,
    }
    cmd[args.command](args)


if __name__ == "__main__":
    main()
