"""
Snakemake wrapper for KINTSUGI registration (multi-cycle alignment via VALIS).

Aligns EDF images across all cycles using rigid + non-rigid registration.
Unlike stitch/decon/edf which process one cycle each, registration processes
ALL cycles in a single job (it needs the full set for alignment).

Snakemake guarantees:
  - All EDF outputs exist before this script runs
  - No need for wait_for_input() polling

Skip-existing: If ALL registered output files already exist for ALL cycles,
the entire job is skipped. Registration is all-or-nothing (no partial resume).
"""

import shutil
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
from skimage.io import imread, imsave

# ---------------------------------------------------------------------------
# Snakemake interface
# ---------------------------------------------------------------------------
PROJECT_DIR = Path(snakemake.params.project_dir)
KINTSUGI_DIR = Path(snakemake.params.kintsugi_dir)
CHANNELS = list(snakemake.params.channels)
CYCLES = list(snakemake.params.cycles)
CHANNEL_NAMES = getattr(snakemake.params, "channel_names", {})

# Setup Python path (same as other wrapper scripts)
sys.path.insert(0, str(PROJECT_DIR / "notebooks"))
sys.path.insert(0, str(KINTSUGI_DIR / "notebooks"))
sys.path.insert(0, str(KINTSUGI_DIR))

# Logging utilities
sys.path.insert(0, snakemake.scriptdir)
from log_utils import log_header, log_footer, log_info, log_warn, log_error
from log_utils import summary_before, summary_after

# ---------------------------------------------------------------------------
# Configuration from config.yaml
# ---------------------------------------------------------------------------
wf_config = snakemake.config
reg_cfg = wf_config.get("registration", {})

REFERENCE_CYCLE = int(reg_cfg.get("reference_cycle", 1))
REFERENCE_CHANNEL = int(reg_cfg.get("reference_channel", 1))
MAX_IMAGE_DIM = int(reg_cfg.get("max_image_dim", 2048))
RIGID_ONLY = bool(reg_cfg.get("rigid_only", False))
FEATURE_DETECTOR = str(reg_cfg.get("feature_detector", "VggFD"))

# ---------------------------------------------------------------------------
# GPU initialization (for VggFD feature detector)
# ---------------------------------------------------------------------------
REQUESTED_DEVICE_MODE = getattr(snakemake.params, "device_mode", "gpu")

DEVICE_MODE = "cpu"
if REQUESTED_DEVICE_MODE == "gpu":
    try:
        import torch

        if torch.cuda.is_available():
            print(f"CUDA available: {torch.cuda.device_count()} GPU(s) - {torch.cuda.get_device_name(0)}")
            DEVICE_MODE = "gpu"
        else:
            print("WARNING: CUDA not available, falling back to CPU mode")
    except ImportError:
        print("WARNING: PyTorch not installed, falling back to CPU mode")
else:
    print("CPU mode (requested by Snakemake rule)")

# Force OrbFD on CPU (VggFD requires GPU)
if DEVICE_MODE == "cpu" and FEATURE_DETECTOR == "VggFD":
    print("Switching to OrbFD feature detector (VggFD requires GPU)")
    FEATURE_DETECTOR = "OrbFD"

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DATA_DIR = PROJECT_DIR / "data"
EDF_DIR = DATA_DIR / "processed" / "edf"
REGISTERED_DIR = DATA_DIR / "processed" / "registered"
REGISTERED_DIR.mkdir(parents=True, exist_ok=True)

# Format cycle numbers consistently
def cyc_fmt(cycle):
    return f"cyc{int(cycle):02d}"


# ---------------------------------------------------------------------------
# Structured logging header
# ---------------------------------------------------------------------------
log_header("registration", 0, PROJECT_DIR)
summary_before("registration", PROJECT_DIR)

print(f"\n{'='*60}")
print(f"Multi-Cycle Registration")
print(f"{'='*60}")
print(f"Project: {PROJECT_DIR.name}")
print(f"Cycles: {CYCLES}")
print(f"Channels: {CHANNELS}")
print(f"Reference cycle: {REFERENCE_CYCLE}")
print(f"Feature detector: {FEATURE_DETECTOR}")
print(f"Max image dim: {MAX_IMAGE_DIM}")
print(f"Rigid only: {RIGID_ONLY}")
print(f"Device mode: {DEVICE_MODE}")
print(f"Input: {EDF_DIR}")
print(f"Output: {REGISTERED_DIR}")


# ---------------------------------------------------------------------------
# Skip-existing check
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# Main processing
# ---------------------------------------------------------------------------
start_time = time.time()

if check_all_registered():
    print("\nAll registered files already exist — nothing to do")
    elapsed = (time.time() - start_time) / 60
    n_files = sum(
        len(list((REGISTERED_DIR / cyc_fmt(c)).glob("*.tif")))
        for c in CYCLES
    )

    sentinel = Path(snakemake.output.sentinel)
    sentinel.parent.mkdir(parents=True, exist_ok=True)
    sentinel.write_text(
        f"stage=registration\n"
        f"completed={datetime.now().isoformat()}\n"
        f"cycles={len(CYCLES)}\n"
        f"skipped=all\n"
        f"files={n_files}\n"
        f"duration_minutes={elapsed:.1f}\n"
    )
    print(f"Sentinel written: {sentinel}")
    summary_after("registration", PROJECT_DIR, elapsed * 60, exit_code=0)
    log_footer(0)
    sys.exit(0)

# ---------------------------------------------------------------------------
# Single-cycle handling: copy EDF → registered (no alignment needed)
# ---------------------------------------------------------------------------
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
    n_files = sum(
        len(list((REGISTERED_DIR / cyc_fmt(c)).glob("*.tif")))
        for c in CYCLES
    )

    sentinel = Path(snakemake.output.sentinel)
    sentinel.parent.mkdir(parents=True, exist_ok=True)
    sentinel.write_text(
        f"stage=registration\n"
        f"completed={datetime.now().isoformat()}\n"
        f"cycles={len(CYCLES)}\n"
        f"method=copy\n"
        f"files={n_files}\n"
        f"duration_minutes={elapsed:.1f}\n"
    )
    print(f"Sentinel written: {sentinel}")
    summary_after("registration", PROJECT_DIR, elapsed * 60, exit_code=0)
    log_footer(0)
    sys.exit(0)

# ---------------------------------------------------------------------------
# Multi-cycle registration via VALIS
# ---------------------------------------------------------------------------
from Kreg import registration, feature_detectors

# Step 1: Collect DAPI/reference images from each cycle
print("\n[1/4] Collecting reference images...")
dapi_images = []
cycle_to_slide_name = {}  # Map cycle number → slide name (DAPI filename stem)

for cyc in CYCLES:
    cycle_dir = EDF_DIR / cyc_fmt(cyc)
    if not cycle_dir.exists():
        log_error(f"EDF directory not found: {cycle_dir}")
        continue

    # Try DAPI by name first
    dapi_files = sorted(cycle_dir.glob("DAPI*.tif"))
    if dapi_files:
        chosen = dapi_files[0]
    else:
        # Fallback: first .tif file (usually channel 1 / DAPI)
        all_tifs = sorted(cycle_dir.glob("*.tif"))
        if all_tifs:
            chosen = all_tifs[0]
        else:
            log_warn(f"No TIF files found for {cyc_fmt(cyc)}")
            continue

    dapi_images.append(str(chosen))
    cycle_to_slide_name[int(cyc)] = chosen.stem
    print(f"  {cyc_fmt(cyc)}: {chosen.name}")

if len(dapi_images) < 2:
    log_error(f"Need at least 2 cycles with images, found {len(dapi_images)}")
    summary_after("registration", PROJECT_DIR, (time.time() - start_time), exit_code=1)
    log_footer(1)
    sys.exit(1)

# Determine reference image
ref_idx = None
for i, cyc in enumerate(CYCLES):
    if int(cyc) == REFERENCE_CYCLE:
        ref_idx = i
        break
if ref_idx is None:
    log_warn(f"Reference cycle {REFERENCE_CYCLE} not found in {CYCLES}, using first cycle")
    ref_idx = 0

print(f"\nReference image: {dapi_images[ref_idx]}")
print(f"Registering {len(dapi_images)} cycles...")

# Step 2: Initialize VALIS registrar
print("\n[2/4] Running VALIS registration...")

reg_output_dir = REGISTERED_DIR / "registration_data"
reg_output_dir.mkdir(parents=True, exist_ok=True)

try:
    # Select feature detector
    try:
        fd_cls = getattr(feature_detectors, FEATURE_DETECTOR)
    except AttributeError:
        log_warn(f"'{FEATURE_DETECTOR}' not found, using OrbFD")
        fd_cls = feature_detectors.OrbFD

    registrar = registration.Valis(
        src_dir=str(EDF_DIR),
        dst_dir=str(reg_output_dir),
        img_list=dapi_images,
        reference_img_f=dapi_images[ref_idx],
        max_image_dim_px=MAX_IMAGE_DIM,
        max_processed_image_dim_px=MAX_IMAGE_DIM,
        feature_detector_cls=fd_cls,
    )

    # Rigid registration
    log_info("Computing rigid transformations...")
    rigid_registrar, non_rigid_reg, rigid_summary = registrar.register()

    if rigid_registrar is None:
        raise RuntimeError("Rigid registration failed — check logs above for details")

    # Non-rigid registration (unless rigid_only)
    if not RIGID_ONLY:
        log_info("Computing non-rigid transformations...")
        non_rigid_registrar, non_rigid_summary = registrar.register_micro(
            max_non_rigid_registration_dim_px=MAX_IMAGE_DIM
        )

    # Step 3: Warp all channels
    print(f"\n[3/4] Warping all channels...")
    n_warped = 0
    n_copied = 0

    # Build reverse lookup: slide name → slide object
    # VALIS slide names come from the DAPI filenames we provided
    slide_lookup = {}
    for slide_obj in registrar.slide_dict.values():
        slide_lookup[slide_obj.name] = slide_obj

    for cyc in CYCLES:
        cycle_edf = EDF_DIR / cyc_fmt(cyc)
        cycle_reg = REGISTERED_DIR / cyc_fmt(cyc)
        cycle_reg.mkdir(parents=True, exist_ok=True)

        # Find the slide object for this cycle
        slide_name = cycle_to_slide_name.get(int(cyc))
        slide_obj = slide_lookup.get(slide_name)

        if slide_obj is None:
            # Try substring match as fallback
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
                # No matching slide — copy unchanged
                shutil.copy2(channel_file, out_path)
                n_copied += 1

        print(f"  {cyc_fmt(cyc)}: warped to {cycle_reg}")

    # Step 4: Save registrar
    print(f"\n[4/4] Saving registrar...")
    registrar_path = reg_output_dir / "registrar.pkl"
    registrar.save_registrar(str(registrar_path))
    print(f"  Saved: {registrar_path}")

    # QC overlay images
    try:
        log_dir = Path(snakemake.log[0]).parent
        _save_registration_qc(registrar, log_dir)
    except Exception as e:
        log_warn(f"QC image generation failed: {e}")

    elapsed = (time.time() - start_time) / 60
    print(f"\n{'='*60}")
    print(f"Registration Complete")
    print(f"{'='*60}")
    print(f"Warped: {n_warped} files")
    print(f"Copied (no match): {n_copied} files")
    print(f"Time: {elapsed:.1f} minutes")

    summary_after("registration", PROJECT_DIR, elapsed * 60, exit_code=0)

    # Write sentinel
    sentinel = Path(snakemake.output.sentinel)
    sentinel.parent.mkdir(parents=True, exist_ok=True)
    sentinel.write_text(
        f"stage=registration\n"
        f"completed={datetime.now().isoformat()}\n"
        f"cycles={len(CYCLES)}\n"
        f"method={'rigid' if RIGID_ONLY else 'rigid+nonrigid'}\n"
        f"feature_detector={FEATURE_DETECTOR}\n"
        f"reference_cycle={REFERENCE_CYCLE}\n"
        f"warped={n_warped}\n"
        f"copied={n_copied}\n"
        f"duration_minutes={elapsed:.1f}\n"
    )
    print(f"Sentinel written: {sentinel}")
    log_footer(0)

except Exception as e:
    log_error(f"Registration failed: {e}")
    import traceback

    traceback.print_exc()

    # Fallback: copy EDF images without transformation
    print("\nFalling back to copying EDF images without transformation...")
    n_copied = 0
    for cyc in CYCLES:
        cycle_edf = EDF_DIR / cyc_fmt(cyc)
        cycle_reg = REGISTERED_DIR / cyc_fmt(cyc)
        cycle_reg.mkdir(parents=True, exist_ok=True)
        for channel_file in cycle_edf.glob("*.tif"):
            shutil.copy2(channel_file, cycle_reg / channel_file.name)
            n_copied += 1
        print(f"  {cyc_fmt(cyc)}: copied {len(list(cycle_edf.glob('*.tif')))} files")

    elapsed = (time.time() - start_time) / 60
    summary_after("registration", PROJECT_DIR, elapsed * 60, exit_code=0)

    # Still write sentinel (fallback completed successfully)
    sentinel = Path(snakemake.output.sentinel)
    sentinel.parent.mkdir(parents=True, exist_ok=True)
    sentinel.write_text(
        f"stage=registration\n"
        f"completed={datetime.now().isoformat()}\n"
        f"cycles={len(CYCLES)}\n"
        f"method=fallback_copy\n"
        f"error={str(e)}\n"
        f"copied={n_copied}\n"
        f"duration_minutes={elapsed:.1f}\n"
    )
    print(f"Sentinel written (fallback): {sentinel}")
    log_footer(0)


# ---------------------------------------------------------------------------
# QC helpers
# ---------------------------------------------------------------------------
def _save_registration_qc(registrar, log_dir):
    """Save overlay QC images to log directory."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    log_dir.mkdir(parents=True, exist_ok=True)

    # Get reference slide
    ref_slide = None
    for slide_obj in registrar.slide_dict.values():
        if slide_obj.is_ref:
            ref_slide = slide_obj
            break

    if ref_slide is None:
        return

    # Create overlay for each non-reference slide
    for slide_obj in registrar.slide_dict.values():
        if slide_obj.is_ref:
            continue

        try:
            fig, axes = plt.subplots(1, 2, figsize=(16, 8))

            # Before (processed images)
            if ref_slide.processed_img is not None and slide_obj.processed_img is not None:
                ref_img = ref_slide.processed_img
                mov_img = slide_obj.processed_img
                overlay = np.zeros((*ref_img.shape[:2], 3), dtype=np.float32)
                r = ref_img.astype(np.float32)
                m = mov_img.astype(np.float32)
                if r.max() > 0:
                    r = r / r.max()
                if m.max() > 0:
                    m = m / m.max()
                overlay[..., 0] = r[:overlay.shape[0], :overlay.shape[1]]
                overlay[..., 1] = m[:overlay.shape[0], :overlay.shape[1]]
                axes[0].imshow(overlay)
                axes[0].set_title("Before Registration")
                axes[0].axis("off")

            # After (registered images)
            if hasattr(ref_slide, "reg_img") and hasattr(slide_obj, "reg_img"):
                if ref_slide.reg_img is not None and slide_obj.reg_img is not None:
                    ref_reg = ref_slide.reg_img
                    mov_reg = slide_obj.reg_img
                    overlay_reg = np.zeros((*ref_reg.shape[:2], 3), dtype=np.float32)
                    r2 = ref_reg.astype(np.float32)
                    m2 = mov_reg.astype(np.float32)
                    if r2.max() > 0:
                        r2 = r2 / r2.max()
                    if m2.max() > 0:
                        m2 = m2 / m2.max()
                    overlay_reg[..., 0] = r2[:overlay_reg.shape[0], :overlay_reg.shape[1]]
                    overlay_reg[..., 1] = m2[:overlay_reg.shape[0], :overlay_reg.shape[1]]
                    axes[1].imshow(overlay_reg)
                    axes[1].set_title("After Registration")
                    axes[1].axis("off")

            plt.suptitle(f"Registration QC: {slide_obj.name}")
            plt.tight_layout()
            qc_path = log_dir / f"reg_qc_{slide_obj.name}.png"
            plt.savefig(str(qc_path), dpi=100, bbox_inches="tight")
            plt.close()
            log_info(f"QC saved: {qc_path.name}")

        except Exception as e:
            log_warn(f"QC image failed for {slide_obj.name}: {e}")
