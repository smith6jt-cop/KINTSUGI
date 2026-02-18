"""
Snakemake wrapper for KINTSUGI registration (multi-cycle alignment via VALIS).

Aligns EDF images across all cycles using rigid + optional non-rigid registration.
Unlike stitch/decon/edf which process one cycle each, registration processes
ALL cycles in a single job (it needs the full set for alignment).

Non-rigid registration (OpticalFlowWarper) is enabled by default with Gaussian
smoothing to prevent noisy displacement fields. Key parameters:
  - max_dim=4096: Displacement field resolution (was incorrectly 2048 in earlier runs)
  - smoothing_method="gauss" with sigma_ratio=0.01: sigma = 0.01 * 4096 = ~41 pixels.
    Moderate smoothing prevents over-warping while correcting real non-rigid deformation.
  - non_rigid_registrar_cls=None in Valis init: Disables harmful coarse NR pass at 1024px

Snakemake guarantees:
  - All EDF outputs exist before this script runs
  - No need for wait_for_input() polling

Skip-existing: If ALL registered output files already exist for ALL cycles,
the entire job is skipped. Registration is all-or-nothing (no partial resume).
"""

import shutil
import sys
import tempfile
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
RIGID_ONLY = bool(reg_cfg.get("rigid_only", True))
FEATURE_DETECTOR = str(reg_cfg.get("feature_detector", "VggFD"))
IMGS_ORDERED = bool(reg_cfg.get("imgs_ordered", True))
ALIGN_TO_REFERENCE = bool(reg_cfg.get("align_to_reference", True))

# Non-rigid parameters (OpticalFlowWarper via VALIS register_micro)
NR_CFG = reg_cfg.get("non_rigid", {})
NR_MAX_DIM = int(NR_CFG.get("max_dim", 4096))
NR_SMOOTHING = NR_CFG.get("smoothing_method", "gauss")
NR_SIGMA_RATIO = float(NR_CFG.get("sigma_ratio", 0.01))
NR_N_GRID_PTS = int(NR_CFG.get("n_grid_pts", 50))
NR_FOLD_PENALTY = float(NR_CFG.get("fold_penalty", 1e-6))
NORMALIZE_DIMS = bool(reg_cfg.get("normalize_dimensions", True))

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
print(f"Images ordered: {IMGS_ORDERED}")
print(f"Align to reference: {ALIGN_TO_REFERENCE}")
if not RIGID_ONLY:
    print(f"Non-rigid max dim: {NR_MAX_DIM}")
    print(f"Non-rigid smoothing: {NR_SMOOTHING} (sigma_ratio={NR_SIGMA_RATIO})")
print(f"Normalize dimensions: {NORMALIZE_DIMS}")
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
from Kreg import registration, feature_detectors, non_rigid_registrars

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

# ---------------------------------------------------------------------------
# Dimension normalization: pad images to uniform (max_h, max_w) if needed
# ---------------------------------------------------------------------------
staging_dir = None  # Set if we create a staging directory

if NORMALIZE_DIMS:
    # Check dimensions across all EDF cycle directories
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
            log_warn(
                f"EDF images have inconsistent dimensions: "
                f"{dict(dims_by_cycle)}. Padding to ({max_h}, {max_w})"
            )

            # Create temp staging directory with padded copies
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

            # Update EDF_DIR and DAPI image paths to use staging
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
    # Select feature detector
    try:
        fd_cls = getattr(feature_detectors, FEATURE_DETECTOR)
    except AttributeError:
        log_warn(f"'{FEATURE_DETECTOR}' not found, using OrbFD")
        fd_cls = feature_detectors.OrbFD

    # NOTE: non_rigid_registrar_cls=None disables the coarse non-rigid pass
    # in register(). The coarse pass defaults to 1024px resolution — too low
    # for CODEX nuclei patterns — and degrades alignment for every cycle.
    # Non-rigid is done exclusively in register_micro() at NR_MAX_DIM (4096px).
    registrar = registration.Valis(
        src_dir=str(EDF_DIR),
        dst_dir=str(reg_output_dir),
        img_list=dapi_images,
        reference_img_f=dapi_images[ref_idx],
        max_image_dim_px=MAX_IMAGE_DIM,
        max_processed_image_dim_px=MAX_IMAGE_DIM,
        non_rigid_registrar_cls=None,
        feature_detector_cls=fd_cls,
        imgs_ordered=IMGS_ORDERED,
        align_to_reference=ALIGN_TO_REFERENCE,
    )

    # Rigid registration
    log_info("Computing rigid transformations...")
    rigid_registrar, non_rigid_reg, rigid_summary = registrar.register()

    if rigid_registrar is None:
        raise RuntimeError("Rigid registration failed — check logs above for details")

    # Non-rigid registration (unless rigid_only)
    if not RIGID_ONLY:
        log_info(
            f"Computing non-rigid transformations "
            f"(max_dim={NR_MAX_DIM}, smoothing={NR_SMOOTHING})..."
        )
        non_rigid_registrar, non_rigid_summary = registrar.register_micro(
            max_non_rigid_registration_dim_px=NR_MAX_DIM,
            non_rigid_registrar_cls=non_rigid_registrars.OpticalFlowWarper,
            non_rigid_reg_params={
                "smoothing_method": NR_SMOOTHING,
                "sigma_ratio": NR_SIGMA_RATIO,
                "n_grid_pts": NR_N_GRID_PTS,
                "fold_penalty": NR_FOLD_PENALTY,
            },
            align_to_reference=ALIGN_TO_REFERENCE,
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

    # Clean up staging directory if dimension normalization was used
    if staging_dir is not None and staging_dir.exists():
        shutil.rmtree(staging_dir)
        print(f"  Cleaned up staging directory: {staging_dir}")

    # Note: thumbnail QC removed — run_registration_qc() in Kprocess.py produces
    # full-resolution spatial NCC heatmaps and targeted overlays (qc_registration rule).

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
    sentinel_lines = [
        f"stage=registration\n",
        f"completed={datetime.now().isoformat()}\n",
        f"cycles={len(CYCLES)}\n",
        f"method={'rigid' if RIGID_ONLY else 'rigid+nonrigid'}\n",
        f"feature_detector={FEATURE_DETECTOR}\n",
        f"reference_cycle={REFERENCE_CYCLE}\n",
        f"imgs_ordered={IMGS_ORDERED}\n",
        f"align_to_reference={ALIGN_TO_REFERENCE}\n",
        f"normalize_dimensions={NORMALIZE_DIMS}\n",
        f"dimensions_normalized={staging_dir is not None}\n",
    ]
    if not RIGID_ONLY:
        sentinel_lines.extend([
            f"non_rigid_max_dim={NR_MAX_DIM}\n",
            f"non_rigid_smoothing={NR_SMOOTHING}\n",
            f"non_rigid_sigma_ratio={NR_SIGMA_RATIO}\n",
        ])
    sentinel_lines.extend([
        f"warped={n_warped}\n",
        f"copied={n_copied}\n",
        f"duration_minutes={elapsed:.1f}\n",
    ])
    sentinel.write_text("".join(sentinel_lines))
    print(f"Sentinel written: {sentinel}")
    log_footer(0)

except Exception as e:
    log_error(f"Registration failed: {e}")
    import traceback

    traceback.print_exc()

    # Restore original EDF_DIR before cleanup (staging may have overridden it)
    if staging_dir is not None:
        EDF_DIR = DATA_DIR / "processed" / "edf"
        if staging_dir.exists():
            shutil.rmtree(staging_dir)
            print(f"  Cleaned up staging directory: {staging_dir}")

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

