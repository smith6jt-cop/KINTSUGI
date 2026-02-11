#!/bin/bash
# =============================================================================
# KINTSUGI Extended Depth of Focus (EDF) Job
# Multi-GPU EDF processing with quality gate and logging
# =============================================================================

# Source utilities (use KINTSUGI_DIR from environment, not relative path)
KINTSUGI_SLURM="${KINTSUGI_DIR}/slurm"
if [ -f "${KINTSUGI_SLURM}/utils.sh" ]; then
    source "${KINTSUGI_SLURM}/utils.sh"
else
    # Fallback: define minimal stubs if utils.sh not found
    log_info() { echo "[INFO] $*"; }
    log_error() { echo "[ERROR] $*" >&2; }
    init_logging() { :; }
    summary_before() { :; }
    summary_after() { :; }
    log_footer() { :; }
    generate_run_id() { date '+%Y%m%d_%H%M%S'; }
fi

# Initialize logging
init_logging "edf" "${RUN_ID:-$(generate_run_id)}"

# Record start time
START_TIME=$(date +%s)

# Generate before summary
summary_before "edf"

echo ""
log_info "Starting EDF processing..."

module load conda
conda activate ${CONDA_ENV:-kintsugi}

export PYTHONPATH="${KINTSUGI_DIR}:${PROJECT_DIR}/notebooks:${PYTHONPATH}"

python << 'PYTHON_SCRIPT'
import sys
import os
import json
import time
import gc
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from itertools import cycle as iter_cycle

PROJECT_DIR = Path(os.environ['PROJECT_DIR'])
KINTSUGI_DIR = Path(os.environ['KINTSUGI_DIR'])
QC_DIR = Path(os.environ.get('QC_DIR', PROJECT_DIR / 'slurm' / 'qc' / 'edf'))
sys.path.insert(0, str(PROJECT_DIR / 'notebooks'))
sys.path.insert(0, str(KINTSUGI_DIR))

from kintsugi.edf import EDFProcessor
from kintsugi.gpu import cleanup_gpu_memory
import numpy as np
import errno
import shutil

# =============================================================================
# METADATA LOADING
# =============================================================================

def load_experiment_config(project_dir):
    """Load experiment configuration from /meta/experiment.json."""
    config_path = project_dir / 'meta' / 'experiment.json'
    if config_path.exists():
        try:
            with open(config_path) as f:
                return json.load(f)
        except (json.JSONDecodeError, IOError) as e:
            print(f"Warning: Could not load experiment.json: {e}")
    return {}

# Load metadata
experiment_config = load_experiment_config(PROJECT_DIR)

print("=" * 60)
print("Metadata Loading")
print("=" * 60)
if experiment_config:
    print(f"  Loaded experiment.json")
else:
    print(f"  experiment.json not found, using environment defaults")

# Check device mode from environment (set by submit.sh based on account type)
DEVICE_MODE = os.environ.get('KINTSUGI_DEVICE_MODE', 'gpu')
print(f"  Device mode: {DEVICE_MODE}")

# Initialize GPU if in GPU mode
if DEVICE_MODE != 'cpu':
    try:
        import cupy as cp
        cp.cuda.Device(0).use()
        _ = cp.zeros(1)  # Test GPU access
        n_gpus = cp.cuda.runtime.getDeviceCount()
        GPU_IDS = list(range(n_gpus))
        print(f"  CUDA initialized: {n_gpus} GPU(s) available")
    except Exception as e:
        print(f"WARNING: CUDA initialization failed: {e}")
        print("Falling back to CPU mode")
        DEVICE_MODE = 'cpu'
        n_gpus = 0
        GPU_IDS = []
else:
    print("  Running in CPU mode (CPU-only burst account)")
    n_gpus = 0
    GPU_IDS = []

# Get configuration: experiment.json -> environment variable -> default
# CYCLE can be passed via --export or derived from SLURM_ARRAY_TASK_ID
CYCLE = int(os.environ.get('CYCLE', os.environ.get('SLURM_ARRAY_TASK_ID', 1)))
START_CHANNEL = int(os.environ.get('START_CHANNEL', 1))
END_CHANNEL = int(os.environ.get('END_CHANNEL',
    experiment_config.get('channels_per_cycle', 4)))
OUTPUT_FORMAT = os.environ.get('OUTPUT_FORMAT', 'zarr')

print(f"  Channels: {START_CHANNEL}-{END_CHANNEL}")

# Load channel names from CHANNELNAMES.txt (matches notebook process_edf_tiff naming)
from Kio import load_channel_names
channels_per_cycle = experiment_config.get('channels_per_cycle', END_CHANNEL)
channel_name_dict = load_channel_names(
    PROJECT_DIR / 'meta',
    channels_per_cycle=channels_per_cycle
)
if channel_name_dict:
    cycle_names = channel_name_dict.get(CYCLE, [])
    print(f"  Channel names (cycle {CYCLE}): {cycle_names}")
else:
    cycle_names = []
    print(f"  WARNING: CHANNELNAMES.txt not found, using CH# fallback naming")

def get_channel_name(ch):
    """Get marker name for channel, with CH# fallback."""
    idx = ch - 1
    if cycle_names and idx < len(cycle_names):
        return cycle_names[idx]
    return f"CH{ch}"

print("=" * 60)

DATA_DIR = PROJECT_DIR / 'data'
DECON_DIR = DATA_DIR / 'processed' / 'deconvolved'
EDF_DIR = DATA_DIR / 'processed' / 'edf'
EDF_DIR.mkdir(parents=True, exist_ok=True)
QC_DIR.mkdir(parents=True, exist_ok=True)

# =============================================================================
# INPUT VALIDATION - Wait for deconvolution to complete
# =============================================================================

def has_data_files(input_dir: Path) -> bool:
    """Check if directory has actual data (TIF files in channel subdirs)."""
    for ch_dir in input_dir.glob("CH*"):
        if ch_dir.is_dir() and any(ch_dir.glob("*.tif")):
            return True
    return False

def wait_for_input(input_dir: Path, max_wait: int = 3600, poll_interval: int = 30) -> bool:
    """
    Wait for input to be ready with polling.

    Checks for a .complete marker file that indicates the upstream job
    finished successfully and all files are written to disk/NFS.
    Also accepts pre-existing data without a marker (e.g., deconvolved before
    markers were implemented).

    Args:
        input_dir: Directory to check for .complete marker
        max_wait: Maximum seconds to wait (default 60 min)
        poll_interval: Seconds between checks (default 30s)

    Returns:
        True if ready, False if timeout
    """
    marker = input_dir / ".complete"
    elapsed = 0
    print(f"Checking input readiness: {input_dir}")

    # Quick check: marker exists
    if marker.exists():
        print(f"Input ready immediately (marker found)")
        return True

    # Quick check: data already present (pre-existing, no marker)
    if has_data_files(input_dir):
        print(f"Input ready immediately (data files found, no marker)")
        return True

    print(f"Waiting for completion marker (timeout: {max_wait}s)...")
    while elapsed < max_wait:
        time.sleep(poll_interval)
        elapsed += poll_interval

        # NFS cache invalidation - list directory to refresh metadata
        try:
            list(input_dir.iterdir())
        except Exception:
            pass

        if marker.exists():
            print(f"Input ready after {elapsed}s")
            return True

        # Also check for data appearing without marker
        if has_data_files(input_dir):
            print(f"Input ready after {elapsed}s (data files found)")
            return True

        if elapsed % 300 == 0:
            print(f"Still waiting ({elapsed}/{max_wait}s)...")

    return False

# Validate deconvolution input is ready (wait up to 60 min for NFS delays / slow upstream)
input_dir = DECON_DIR / f"cyc{CYCLE:02d}"
if not input_dir.exists():
    print(f"ERROR: Deconvolution output directory does not exist: {input_dir}")
    print("The deconvolution job may not have started or failed early.")
    sys.exit(1)

if not wait_for_input(input_dir, max_wait=3600, poll_interval=30):
    print(f"ERROR: Deconvolution incomplete for cycle {CYCLE}")
    print(f"Missing completion marker: {input_dir / '.complete'}")
    print(f"Check deconvolution job status: squeue -u $USER")
    print(f"Check deconvolution logs for errors.")
    sys.exit(1)

print(f"Deconvolution input validated for cycle {CYCLE}")

# =============================================================================
# CYCLE CLAIMING (prevents GPU/CPU duplicate work)
# =============================================================================
def claim_cycle(lock_base, cycle):
    """Atomically claim a cycle using mkdir (atomic on POSIX/NFS)."""
    lock_dir = lock_base / f".kintsugi_lock_cyc{cycle:02d}"
    my_job_id = os.environ.get('SLURM_JOB_ID', 'unknown')
    try:
        os.mkdir(str(lock_dir))
        with open(str(lock_dir / "owner.txt"), 'w') as f:
            f.write(f"job_id={my_job_id}\n")
            f.write(f"device_mode={os.environ.get('KINTSUGI_DEVICE_MODE', 'unknown')}\n")
            f.write(f"hostname={os.environ.get('HOSTNAME', 'unknown')}\n")
        return True
    except FileExistsError:
        return False
    except OSError as e:
        if e.errno == errno.EEXIST:
            return False
        raise

def check_lock_owner(lock_base, cycle):
    """Check if existing lock belongs to this job (requeue case)."""
    lock_dir = lock_base / f".kintsugi_lock_cyc{cycle:02d}"
    owner_file = lock_dir / "owner.txt"
    my_job_id = os.environ.get('SLURM_JOB_ID', 'unknown')
    if owner_file.exists():
        try:
            content = owner_file.read_text()
            for line in content.strip().split('\n'):
                if line.startswith('job_id=') and line.split('=', 1)[1] == my_job_id:
                    return True
        except Exception:
            pass
    return False

def release_lock(lock_base, cycle):
    """Remove lock after successful processing."""
    lock_dir = lock_base / f".kintsugi_lock_cyc{cycle:02d}"
    shutil.rmtree(str(lock_dir), ignore_errors=True)

def cleanup_stale_locks(lock_base, max_age_hours=12):
    """Remove locks older than max_age_hours (handles dead jobs)."""
    import time as _time
    for lock_dir in lock_base.glob(".kintsugi_lock_cyc*"):
        owner_file = lock_dir / "owner.txt"
        if owner_file.exists():
            try:
                age_hours = (_time.time() - owner_file.stat().st_mtime) / 3600
                if age_hours > max_age_hours:
                    shutil.rmtree(str(lock_dir), ignore_errors=True)
                    print(f"  Cleaned stale lock: {lock_dir.name} (age: {age_hours:.1f}h)")
            except Exception:
                pass

# =============================================================================
# EDF PARAMETERS - Updated to exclude edge z-planes and prevent distortion
# =============================================================================
EDF_PARAMS = {
    'radius_x': 2,
    'radius_y': 2,
    'sigma': 10.0,
    'tiles': (3, 3),      # Process in 9 tiles to fit in GPU memory
    'backend': 'numpy' if DEVICE_MODE == 'cpu' else 'auto',
    'device': DEVICE_MODE,
    'blend_depth': 0,
    'z_smooth_sigma': 1.0
}

# Quality gate thresholds
QUALITY_GATE = {
    'max_zero_pct': 0.10,    # Exclude slices with >10% zeros
    'max_sat_pct': 0.05,     # Exclude slices with >5% saturation
    'min_valid_slices': 3,   # Need at least 3 valid slices
}

print(f"\n{'='*60}")
print(f"EDF Processing - Cycle {CYCLE}")
print(f"{'='*60}")
print(f"Project: {PROJECT_DIR.name}")
print(f"Channels: {START_CHANNEL}-{END_CHANNEL}")
print(f"Device mode: {DEVICE_MODE}")
if DEVICE_MODE == 'gpu':
    print(f"GPUs: {n_gpus}")
else:
    print("Running on CPU (burst account)")
print(f"Input: {DECON_DIR}")
print(f"Output: {EDF_DIR}")
print(f"QC output: {QC_DIR}")
print(f"EDF params: z exclusion enabled, quality gate enabled")

# =============================================================================
# SKIP-EXISTING CHECK
# =============================================================================
def check_cycle_complete(edf_dir, cycle, start_ch, end_ch):
    """Check if EDF output already exists and is complete."""
    cycle_out = edf_dir / f"cyc{cycle:02d}"
    if not cycle_out.exists():
        return False, "output directory missing"

    for ch in range(start_ch, end_ch + 1):
        ch_name = get_channel_name(ch)
        edf_file = cycle_out / f"{ch_name}.tif"
        if not edf_file.exists():
            return False, f"{ch_name}.tif missing"

    return True, "complete"

# Check if this cycle is already processed
is_complete, status = check_cycle_complete(EDF_DIR, CYCLE, START_CHANNEL, END_CHANNEL)
if is_complete:
    print(f"\n{'='*60}")
    print(f"SKIP: Cycle {CYCLE} EDF already complete ({status})")
    print(f"{'='*60}")
    print(f"Output exists: {EDF_DIR / f'cyc{CYCLE:02d}'}")
    print("To reprocess, delete the output directory first.")
    sys.exit(0)
else:
    print(f"\nCycle {CYCLE} needs EDF processing: {status}")

# Claim this cycle (prevents GPU/CPU pool from duplicating work)
cleanup_stale_locks(EDF_DIR, max_age_hours=12)
if not claim_cycle(EDF_DIR, CYCLE):
    # Check if this is a requeued job reclaiming its own lock
    if check_lock_owner(EDF_DIR, CYCLE):
        print(f"Requeued job reclaiming cycle {CYCLE} lock")
    else:
        print(f"SKIP: Cycle {CYCLE} already claimed by another job")
        sys.exit(0)

channels = list(range(START_CHANNEL, END_CHANNEL + 1))

def load_stack_parallel(tiff_files):
    """Load z-stack from TIFF files with parallel I/O."""
    from concurrent.futures import ThreadPoolExecutor
    from skimage.io import imread

    def load_single(f):
        return imread(str(f))

    with ThreadPoolExecutor(max_workers=4) as executor:
        slices = list(executor.map(load_single, tiff_files))

    return np.stack(slices, axis=0)

def quality_gate_zstack(stack, channel_name):
    """
    Quality gate: detect and exclude corrupted z-slices.

    Returns tuple of (valid_indices, excluded_reasons)
    """
    valid_z_indices = []
    excluded_reasons = []

    for z_idx in range(stack.shape[0]):
        slice_data = stack[z_idx]
        total_pixels = slice_data.size

        pct_zero = np.sum(slice_data == 0) / total_pixels
        pct_sat = np.sum(slice_data >= 64000) / total_pixels

        if pct_zero > QUALITY_GATE['max_zero_pct']:
            excluded_reasons.append(f"z={z_idx+1}: {pct_zero*100:.1f}% zeros")
            continue

        if pct_sat > QUALITY_GATE['max_sat_pct']:
            excluded_reasons.append(f"z={z_idx+1}: {pct_sat*100:.1f}% saturated")
            continue

        valid_z_indices.append(z_idx)

    return valid_z_indices, excluded_reasons

def save_qc_image(data, output_path, title=""):
    """Save EDF result as QC image."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 8))
        vmin, vmax = np.percentile(data, [1, 99])
        im = ax.imshow(data, cmap='gray', vmin=vmin, vmax=vmax)
        ax.set_title(title)
        ax.axis('off')
        plt.colorbar(im, ax=ax, fraction=0.046)
        plt.tight_layout()
        plt.savefig(output_path, dpi=100, bbox_inches='tight')
        plt.close()
        print(f"    QC saved: {output_path.name}")
    except Exception as e:
        print(f"    QC save failed: {e}")

def process_edf(args):
    """Process EDF for one channel with quality gate."""
    ch, gpu_id = args
    ch_name = get_channel_name(ch)
    if DEVICE_MODE == 'gpu':
        print(f"\n  [GPU{gpu_id}] Channel {ch} ({ch_name}) starting...")
    else:
        print(f"\n  [CPU] Channel {ch} ({ch_name}) starting...")

    try:
        if DEVICE_MODE == 'gpu' and n_gpus > 0:
            import cupy as cp
            cp.cuda.Device(gpu_id).use()

        processor = EDFProcessor(backend=EDF_PARAMS['backend'], method='variance')

        input_path = DECON_DIR / f"cyc{CYCLE:02d}" / f"CH{ch}"
        output_path = EDF_DIR / f"cyc{CYCLE:02d}"
        output_path.mkdir(parents=True, exist_ok=True)

        # Find input files
        input_files = sorted(input_path.glob("*.tif"))
        if not input_files:
            print(f"    No input files found for channel {ch}")
            return ch, False, "No input files"

        # Load z-stack with parallel I/O
        print(f"    Loading {len(input_files)} z-slices...")
        stack = load_stack_parallel(input_files)
        n_zplanes = stack.shape[0]
        print(f"    Stack shape: {stack.shape}")

        # ======== QUALITY GATE ========
        valid_z_indices, excluded_reasons = quality_gate_zstack(stack, f"CH{ch}")

        if excluded_reasons:
            print(f"    Quality gate excluded z-slices:")
            for reason in excluded_reasons:
                print(f"      - {reason}")

        if len(valid_z_indices) < QUALITY_GATE['min_valid_slices']:
            msg = f"Only {len(valid_z_indices)} valid z-slices, need at least {QUALITY_GATE['min_valid_slices']}"
            print(f"    ERROR: {msg}")
            return ch, False, msg

        # Filter to valid slices
        stack_filtered = stack[valid_z_indices]
        print(f"    Valid z-slices: {len(valid_z_indices)}/{n_zplanes}")

        # Apply z_start/z_end bounds (exclude first and last of remaining)
        if len(valid_z_indices) > 4:
            z_start = 2  # 1-indexed, skip first valid slice
            z_end = len(valid_z_indices) - 1  # Skip last valid slice
        else:
            z_start = 1
            z_end = len(valid_z_indices)

        print(f"    Using z_start={z_start}, z_end={z_end}")

        # ======== RUN EDF ========
        edf_result = processor.process(
            stack_filtered,
            radius_x=EDF_PARAMS['radius_x'],
            radius_y=EDF_PARAMS['radius_y'],
            sigma=EDF_PARAMS['sigma'],
            z_start=z_start,
            z_end=z_end,
            tiles=EDF_PARAMS['tiles'],
            device=EDF_PARAMS['device'],
            blend_depth=EDF_PARAMS.get('blend_depth', 0),
            z_smooth_sigma=EDF_PARAMS.get('z_smooth_sigma', 0.0)
        )

        # Save result (use marker name to match notebook process_edf_tiff convention)
        from skimage.io import imsave
        output_file = output_path / f"{ch_name}.tif"
        imsave(str(output_file), edf_result.astype(np.uint16), check_contrast=False)

        # Verify output quality
        pct_zero = 100 * np.sum(edf_result == 0) / edf_result.size
        print(f"    Output: {pct_zero:.2f}% zeros (should be <5%)")

        if pct_zero > 5:
            print(f"    WARNING: High zero percentage in output")

        # Generate QC image
        qc_path = QC_DIR / f"cyc{CYCLE:02d}_{ch_name}_edf.png"
        save_qc_image(edf_result, qc_path, f"Cycle {CYCLE} {ch_name} EDF ({pct_zero:.1f}% zeros)")

        if DEVICE_MODE == 'gpu':
            print(f"  [GPU{gpu_id}] Channel {ch} ({ch_name}) complete -> {output_file.name}")
        else:
            print(f"  [CPU] Channel {ch} ({ch_name}) complete -> {output_file.name}")

        # Cleanup
        del stack, stack_filtered, edf_result
        gc.collect()

        return ch, True, None

    except Exception as e:
        if DEVICE_MODE == 'gpu':
            print(f"  [GPU{gpu_id}] Channel {ch} ({ch_name}) FAILED: {e}")
        else:
            print(f"  [CPU] Channel {ch} ({ch_name}) FAILED: {e}")
        import traceback
        traceback.print_exc()
        return ch, False, str(e)
    finally:
        cleanup_gpu_memory()

# Process channels
start_time = time.time()

if DEVICE_MODE == 'cpu':
    # CPU mode: process channels sequentially
    print(f"\n[CPU] Processing {len(channels)} channels sequentially")
    results = [process_edf((ch, 0)) for ch in channels]
elif n_gpus >= 2 and len(channels) >= 2:
    # Multi-GPU mode
    print(f"\n[MULTI-GPU] {len(channels)} channels across {n_gpus} GPUs")
    pairs = list(zip(channels, iter_cycle(GPU_IDS)))
    with ThreadPoolExecutor(max_workers=n_gpus) as executor:
        results = list(executor.map(process_edf, pairs))
else:
    # Single GPU mode
    print(f"\n[GPU] Processing {len(channels)} channels on GPU 0")
    results = [process_edf((ch, 0)) for ch in channels]

successful = sum(1 for _, ok, _ in results if ok)
elapsed = (time.time() - start_time) / 60

print(f"\n{'='*60}")
print(f"EDF Complete")
print(f"{'='*60}")
print(f"Success: {successful}/{len(channels)} channels")
print(f"Time: {elapsed:.1f} minutes")
print(f"Output: {EDF_DIR}")
print(f"QC images: {QC_DIR}")

# Report any failures
failed = [(ch, err) for ch, ok, err in results if not ok]
if failed:
    print(f"\nFailed channels:")
    for ch, err in failed:
        print(f"  {get_channel_name(ch)} (CH{ch}): {err}")
    sys.exit(1)

# Release cycle lock (processing succeeded)
release_lock(EDF_DIR, CYCLE)

# =============================================================================
# WRITE COMPLETION MARKER
# =============================================================================
# This marker signals to downstream jobs (registration) that EDF is complete
# and all output files have been written and flushed to disk/NFS.

from datetime import datetime

output_cycle_dir = EDF_DIR / f"cyc{CYCLE:02d}"
marker_file = output_cycle_dir / ".complete"

try:
    with open(marker_file, 'w') as f:
        f.write(f"stage=edf\n")
        f.write(f"cycle={CYCLE}\n")
        f.write(f"completed={datetime.now().isoformat()}\n")
        f.write(f"channels={START_CHANNEL}-{END_CHANNEL}\n")
        f.write(f"successful={successful}\n")
        f.write(f"job_id={os.environ.get('SLURM_JOB_ID', 'local')}\n")
        f.write(f"array_task_id={os.environ.get('SLURM_ARRAY_TASK_ID', '0')}\n")
        f.write(f"duration_minutes={elapsed:.1f}\n")

    # Force NFS sync to ensure marker is visible to other nodes
    os.sync()
    print(f"Completion marker written: {marker_file}")
except Exception as e:
    print(f"WARNING: Failed to write completion marker: {e}")
    # Don't fail the job if marker write fails - the data is still valid
PYTHON_SCRIPT

EXIT_CODE=$?

echo ""
log_info "EDF processing finished with exit code: ${EXIT_CODE}"

# Generate after summary
summary_after "edf" "${EXIT_CODE}" "${START_TIME}"

# Log footer
log_footer "${EXIT_CODE}"

exit ${EXIT_CODE}
