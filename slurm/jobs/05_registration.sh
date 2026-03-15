#!/bin/bash
# =============================================================================
# KINTSUGI Multi-Cycle Registration Job
# GPU-accelerated VALIS registration (rigid + non-rigid) across all cycles
# =============================================================================
#
# Prerequisites: All EDF outputs complete (data/processed/edf/cyc{NN}/)
# Output: Aligned images in data/processed/registered/cyc{NN}/
#
# Usage:
#   sbatch --export=ALL,PROJECT_DIR=/path/to/project 05_registration.sh
#
# NOTE: Unlike stitch/decon/edf which run per-cycle, registration processes
# ALL cycles in a single job (multi-cycle alignment requires the full set).
# =============================================================================

#SBATCH --job-name=kintsugi_registration
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64gb
#SBATCH --time=04:00:00
#SBATCH --gres=gpu:1
#SBATCH --qos=${QOS:-${SLURM_ACCOUNT}-b}

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
init_logging "registration" "${RUN_ID:-$(generate_run_id)}"

# Record start time
START_TIME=$(date +%s)

# Generate before summary
summary_before "registration"

echo ""
log_info "Starting multi-cycle registration..."

module load conda
conda activate ${CONDA_ENV:-kintsugi}

export PYTHONPATH="${KINTSUGI_DIR}:${KINTSUGI_DIR}/notebooks:${PROJECT_DIR}/notebooks:${PYTHONPATH}"

python << 'PYTHON_SCRIPT'
import sys
import os
import json
import shutil
import time
import tempfile
from pathlib import Path
from datetime import datetime

import numpy as np
from skimage.io import imread, imsave

PROJECT_DIR = Path(os.environ['PROJECT_DIR'])
KINTSUGI_DIR = Path(os.environ['KINTSUGI_DIR'])
sys.path.insert(0, str(PROJECT_DIR / 'notebooks'))
sys.path.insert(0, str(KINTSUGI_DIR / 'notebooks'))
sys.path.insert(0, str(KINTSUGI_DIR))

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

experiment_config = load_experiment_config(PROJECT_DIR)

# =============================================================================
# CONFIGURATION
# =============================================================================

# Registration parameters (from environment or defaults)
REFERENCE_CYCLE = int(os.environ.get('REG_REFERENCE_CYCLE', 1))
MAX_IMAGE_DIM = int(os.environ.get('REG_MAX_IMAGE_DIM', 4096))
RIGID_ONLY = os.environ.get('REG_RIGID_ONLY', 'false').lower() == 'true'
FEATURE_DETECTOR = os.environ.get('REG_FEATURE_DETECTOR', 'VggFD')
NR_MAX_DIM = int(os.environ.get('REG_NR_MAX_DIM', 4096))
NR_SMOOTHING = os.environ.get('REG_NR_SMOOTHING', 'gauss')
NR_SIGMA_RATIO = float(os.environ.get('REG_NR_SIGMA_RATIO', '0.01'))
NORMALIZE_DIMS = os.environ.get('REG_NORMALIZE_DIMS', 'true').lower() == 'true'

# Discover cycles from EDF output
DATA_DIR = PROJECT_DIR / 'data'
EDF_DIR = DATA_DIR / 'processed' / 'edf'
REGISTERED_DIR = DATA_DIR / 'processed' / 'registered'
REGISTERED_DIR.mkdir(parents=True, exist_ok=True)

# Find all completed EDF cycles
CYCLES = sorted([
    int(d.name.replace('cyc', ''))
    for d in EDF_DIR.iterdir()
    if d.is_dir() and d.name.startswith('cyc') and any(d.glob('*.tif'))
])

if not CYCLES:
    print("ERROR: No EDF cycle directories found")
    sys.exit(1)

def cyc_fmt(cycle):
    return f"cyc{int(cycle):02d}"

# =============================================================================
# GPU INITIALIZATION
# =============================================================================
DEVICE_MODE = os.environ.get('KINTSUGI_DEVICE_MODE', 'gpu')

if DEVICE_MODE == 'gpu':
    try:
        import torch
        if torch.cuda.is_available():
            print(f"CUDA available: {torch.cuda.device_count()} GPU(s) - {torch.cuda.get_device_name(0)}")
        else:
            print("WARNING: CUDA not available, falling back to CPU mode")
            DEVICE_MODE = 'cpu'
    except ImportError:
        print("WARNING: PyTorch not installed, falling back to CPU mode")
        DEVICE_MODE = 'cpu'

# Force OrbFD on CPU (VggFD requires GPU)
if DEVICE_MODE == 'cpu' and FEATURE_DETECTOR == 'VggFD':
    print("Switching to OrbFD feature detector (VggFD requires GPU)")
    FEATURE_DETECTOR = 'OrbFD'

print(f"\n{'='*60}")
print(f"Multi-Cycle Registration")
print(f"{'='*60}")
print(f"Project: {PROJECT_DIR.name}")
print(f"Cycles: {CYCLES}")
print(f"Reference cycle: {REFERENCE_CYCLE}")
print(f"Feature detector: {FEATURE_DETECTOR}")
print(f"Max image dim: {MAX_IMAGE_DIM}")
print(f"Rigid only: {RIGID_ONLY}")
if not RIGID_ONLY:
    print(f"Non-rigid max dim: {NR_MAX_DIM}")
    print(f"Non-rigid smoothing: {NR_SMOOTHING} (sigma_ratio={NR_SIGMA_RATIO})")
print(f"Normalize dimensions: {NORMALIZE_DIMS}")
print(f"Device mode: {DEVICE_MODE}")

# =============================================================================
# SKIP-EXISTING CHECK
# =============================================================================

def check_all_registered():
    """Check if ALL registered output files exist for ALL cycles."""
    for cyc in CYCLES:
        cycle_edf = EDF_DIR / cyc_fmt(cyc)
        cycle_reg = REGISTERED_DIR / cyc_fmt(cyc)
        if not cycle_edf.exists():
            return False
        edf_files = list(cycle_edf.glob("*.tif"))
        if not edf_files:
            return False
        for edf_file in edf_files:
            if not (cycle_reg / edf_file.name).exists():
                return False
    return True

start_time = time.time()

if check_all_registered():
    print("\nAll registered files already exist — nothing to do")
    elapsed = (time.time() - start_time) / 60
    print(f"Time: {elapsed:.1f} minutes")
    sys.exit(0)

# =============================================================================
# SINGLE-CYCLE HANDLING
# =============================================================================

if len(CYCLES) < 2:
    print("\nSingle cycle — copying EDF to registered (no alignment needed)")
    for cyc in CYCLES:
        src = EDF_DIR / cyc_fmt(cyc)
        dst = REGISTERED_DIR / cyc_fmt(cyc)
        dst.mkdir(parents=True, exist_ok=True)
        for f in src.glob("*.tif"):
            shutil.copy2(f, dst / f.name)
            print(f"  Copied: {f.name}")

    elapsed = (time.time() - start_time) / 60
    print(f"Time: {elapsed:.1f} minutes")
    sys.exit(0)

# =============================================================================
# MULTI-CYCLE REGISTRATION VIA VALIS (local Kreg)
# =============================================================================

from Kreg import registration, feature_detectors, non_rigid_registrars

# Step 1: Collect DAPI/reference images from each cycle
print("\n[1/4] Collecting reference images...")
dapi_images = []
cycle_to_slide_name = {}

for cyc in CYCLES:
    cycle_dir = EDF_DIR / cyc_fmt(cyc)
    if not cycle_dir.exists():
        print(f"  WARNING: EDF directory not found: {cycle_dir}")
        continue

    dapi_files = sorted(cycle_dir.glob("DAPI*.tif"))
    if dapi_files:
        chosen = dapi_files[0]
    else:
        all_tifs = sorted(cycle_dir.glob("*.tif"))
        if all_tifs:
            chosen = all_tifs[0]
        else:
            print(f"  WARNING: No TIF files found for {cyc_fmt(cyc)}")
            continue

    dapi_images.append(str(chosen))
    cycle_to_slide_name[int(cyc)] = chosen.stem
    print(f"  {cyc_fmt(cyc)}: {chosen.name}")

if len(dapi_images) < 2:
    print(f"ERROR: Need at least 2 cycles with images, found {len(dapi_images)}")
    sys.exit(1)

# Determine reference image
ref_idx = 0
for i, cyc in enumerate(CYCLES):
    if int(cyc) == REFERENCE_CYCLE:
        ref_idx = i
        break

print(f"\nReference image: {dapi_images[ref_idx]}")

# =============================================================================
# DIMENSION NORMALIZATION
# =============================================================================
staging_dir = None

if NORMALIZE_DIMS:
    dims_by_cycle = {}
    for cyc in CYCLES:
        cycle_dir = EDF_DIR / cyc_fmt(cyc)
        if not cycle_dir.exists():
            continue
        first_tif = next(cycle_dir.glob("*.tif"), None)
        if first_tif is not None:
            img = imread(str(first_tif))
            dims_by_cycle[int(cyc)] = img.shape[:2]

    if dims_by_cycle:
        unique_dims = set(dims_by_cycle.values())
        if len(unique_dims) > 1:
            max_h = max(h for h, w in unique_dims)
            max_w = max(w for h, w in unique_dims)
            print(f"  WARNING: Inconsistent dims: {dict(dims_by_cycle)}. Padding to ({max_h}, {max_w})")

            staging_dir = Path(tempfile.mkdtemp(
                prefix="kintsugi_reg_staging_",
                dir=str(DATA_DIR / "processed"),
            ))

            for cyc in CYCLES:
                src_cycle = EDF_DIR / cyc_fmt(cyc)
                dst_cycle = staging_dir / cyc_fmt(cyc)
                dst_cycle.mkdir(parents=True, exist_ok=True)
                if not src_cycle.exists():
                    continue
                h, w = dims_by_cycle.get(int(cyc), (max_h, max_w))
                needs_pad = (h != max_h) or (w != max_w)
                for tif in sorted(src_cycle.glob("*.tif")):
                    if needs_pad:
                        img = imread(str(tif))
                        padded = np.zeros((max_h, max_w), dtype=img.dtype)
                        padded[:img.shape[0], :img.shape[1]] = img
                        imsave(str(dst_cycle / tif.name), padded, check_contrast=False)
                    else:
                        shutil.copy2(tif, dst_cycle / tif.name)

            original_edf_dir = EDF_DIR
            EDF_DIR = staging_dir
            dapi_images = [
                str(staging_dir / Path(p).relative_to(original_edf_dir))
                for p in dapi_images
            ]
            print(f"  Padded images staged to: {staging_dir}")
        else:
            dims = next(iter(unique_dims))
            print(f"  All EDF images have uniform dimensions: {dims}")

# Step 2: Initialize VALIS registrar
print("\n[2/4] Running VALIS registration...")

reg_output_dir = REGISTERED_DIR / "registration_data"
reg_output_dir.mkdir(parents=True, exist_ok=True)

try:
    try:
        fd_cls = getattr(feature_detectors, FEATURE_DETECTOR)
    except AttributeError:
        print(f"  WARNING: '{FEATURE_DETECTOR}' not found, using OrbFD")
        fd_cls = feature_detectors.OrbFD

    # non_rigid_registrar_cls=None disables coarse NR pass at 1024px
    registrar = registration.Valis(
        src_dir=str(EDF_DIR),
        dst_dir=str(reg_output_dir),
        img_list=dapi_images,
        reference_img_f=dapi_images[ref_idx],
        max_image_dim_px=MAX_IMAGE_DIM,
        max_processed_image_dim_px=MAX_IMAGE_DIM,
        non_rigid_registrar_cls=None,
        feature_detector_cls=fd_cls,
        imgs_ordered=True,
        align_to_reference=True,
    )

    print("  Computing rigid transformations...")
    rigid_registrar, non_rigid_reg, rigid_summary = registrar.register()

    if rigid_registrar is None:
        raise RuntimeError("Rigid registration failed — check logs above")

    if not RIGID_ONLY:
        print(f"  Computing non-rigid transformations (max_dim={NR_MAX_DIM}, smoothing={NR_SMOOTHING})...")
        non_rigid_registrar, non_rigid_summary = registrar.register_micro(
            max_non_rigid_registration_dim_px=NR_MAX_DIM,
            non_rigid_registrar_cls=non_rigid_registrars.OpticalFlowWarper,
            non_rigid_reg_params={
                "smoothing_method": NR_SMOOTHING,
                "sigma_ratio": NR_SIGMA_RATIO,
            },
            align_to_reference=True,
        )

    # Step 3: Warp all channels
    print(f"\n[3/4] Warping all channels...")
    n_warped = 0

    slide_lookup = {}
    for slide_obj in registrar.slide_dict.values():
        slide_lookup[slide_obj.name] = slide_obj

    for cyc in CYCLES:
        cycle_edf = EDF_DIR / cyc_fmt(cyc)
        cycle_reg = REGISTERED_DIR / cyc_fmt(cyc)
        cycle_reg.mkdir(parents=True, exist_ok=True)

        slide_name = cycle_to_slide_name.get(int(cyc))
        slide_obj = slide_lookup.get(slide_name)

        if slide_obj is None:
            for sname, sobj in slide_lookup.items():
                if cyc_fmt(cyc) in sname or str(int(cyc)) in sname:
                    slide_obj = sobj
                    break

        for channel_file in sorted(cycle_edf.glob("*.tif")):
            out_path = cycle_reg / channel_file.name
            if slide_obj is not None:
                img = imread(str(channel_file))
                warped = slide_obj.warp_img(img)
                imsave(str(out_path), warped.astype(np.uint16), check_contrast=False)
                n_warped += 1
            else:
                shutil.copy2(channel_file, out_path)

        print(f"  {cyc_fmt(cyc)}: warped to {cycle_reg}")

    # Step 4: Save registrar
    print(f"\n[4/4] Saving registrar...")
    registrar_path = reg_output_dir / "registrar.pkl"
    registrar.save_registrar(str(registrar_path))
    print(f"  Saved: {registrar_path}")

    if staging_dir is not None and staging_dir.exists():
        shutil.rmtree(staging_dir)

    elapsed = (time.time() - start_time) / 60
    print(f"\n{'='*60}")
    print(f"Registration Complete")
    print(f"{'='*60}")
    print(f"Warped: {n_warped} files")
    print(f"Time: {elapsed:.1f} minutes")

    # Write completion marker
    marker_file = REGISTERED_DIR / ".complete"
    with open(marker_file, 'w') as f:
        f.write(f"stage=registration\n")
        f.write(f"completed={datetime.now().isoformat()}\n")
        f.write(f"cycles={len(CYCLES)}\n")
        f.write(f"method={'rigid' if RIGID_ONLY else 'rigid+nonrigid'}\n")
        f.write(f"feature_detector={FEATURE_DETECTOR}\n")
        f.write(f"reference_cycle={REFERENCE_CYCLE}\n")
        f.write(f"warped={n_warped}\n")
        f.write(f"job_id={os.environ.get('SLURM_JOB_ID', 'local')}\n")
        f.write(f"duration_minutes={elapsed:.1f}\n")

    os.sync()
    print(f"Completion marker written: {marker_file}")

except Exception as e:
    print(f"ERROR: Registration failed: {e}")
    import traceback
    traceback.print_exc()

    if staging_dir is not None and staging_dir.exists():
        shutil.rmtree(staging_dir)

    sys.exit(1)
PYTHON_SCRIPT

EXIT_CODE=$?

echo ""
log_info "Registration finished with exit code: ${EXIT_CODE}"

# Generate after summary
summary_after "registration" "${EXIT_CODE}" "${START_TIME}"

# Log footer
log_footer "${EXIT_CODE}"

exit ${EXIT_CODE}
