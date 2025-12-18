"""
Interactive artifact detection and mitigation tuner for Jupyter notebooks.

Provides widget-based interface for reviewing detected artifacts,
selecting mitigation methods, and verifying results.
"""

import json
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

import numpy as np

try:
    import ipywidgets as widgets
    from IPython.display import clear_output, display
    WIDGETS_AVAILABLE = True
except ImportError:
    WIDGETS_AVAILABLE = False

try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


@dataclass
class ArtifactTunerResult:
    """Result from the artifact tuner session."""

    channel_name: str
    decisions: dict = field(default_factory=dict)  # z_index -> (method, mitigated_image)
    skipped_planes: list = field(default_factory=list)  # Planes marked as no artifact
    mitigated_planes: dict = field(default_factory=dict)  # z_index -> mitigated image
    approved: bool = False

    def get_mitigated_image(self, z_index: int) -> np.ndarray | None:
        """Get mitigated image for a z-plane, or None if not mitigated."""
        return self.mitigated_planes.get(z_index)


def artifact_detection_tuner(
    z_stack: np.ndarray,
    channel_name: str = "",
    sensitivity: float = 0.5,
    output_dir: Path | None = None,
) -> ArtifactTunerResult:
    """
    Interactive widget for artifact detection and mitigation.

    Features:
    - Z-plane slider to browse stack
    - FFT visualization showing frequency peaks
    - Mitigation method dropdown
    - Before/After comparison
    - Approve button to confirm decisions

    Parameters
    ----------
    z_stack : np.ndarray
        3D array of shape (nz, height, width)
    channel_name : str
        Name for display/logging
    sensitivity : float
        Detection sensitivity (0-1)
    output_dir : Path, optional
        Directory to save mitigated files and report

    Returns
    -------
    ArtifactTunerResult
        Contains all decisions and mitigated images
    """
    from ..qc.stripe_artifact import (
        detect_stripe_artifacts,
        get_fft_visualization,
        scan_zstack_for_artifacts,
    )
    from ..qc.stripe_mitigation import get_recommended_method, mitigate_artifact

    # Scan for artifacts
    print(f"Scanning {channel_name} for artifacts (sensitivity={sensitivity})...")
    report = scan_zstack_for_artifacts(z_stack, sensitivity=sensitivity, verbose=False)

    if not report.has_artifacts:
        print(f"[OK] No artifacts detected in {channel_name}")
        return ArtifactTunerResult(channel_name=channel_name, approved=True)

    print(f"[WARNING] {report}")

    # Initialize result
    result = ArtifactTunerResult(channel_name=channel_name)

    if not WIDGETS_AVAILABLE or not MATPLOTLIB_AVAILABLE:
        # Fallback to non-interactive mode
        print("\nInteractive widgets not available. Using automatic mitigation.")
        return _auto_mitigate(z_stack, report, result, output_dir)

    # Create interactive widget interface
    nz = z_stack.shape[0]

    # State
    current_decisions = {}
    current_mitigated = {}

    # Widgets
    z_slider = widgets.IntSlider(
        value=report.affected_planes[0] if report.affected_planes else 0,
        min=0,
        max=nz - 1,
        step=1,
        description='Z-plane:',
        continuous_update=False,
    )

    method_dropdown = widgets.Dropdown(
        options=[
            ('Notch Filter (Recommended)', 'notch'),
            ('Directional Blur', 'directional'),
            ('Interpolate from neighbors', 'interpolate'),
            ('Keep as-is (no artifact)', 'keep_as_is'),
        ],
        value='notch',
        description='Method:',
    )

    sensitivity_slider = widgets.FloatSlider(
        value=sensitivity,
        min=0.1,
        max=1.0,
        step=0.1,
        description='Sensitivity:',
        continuous_update=False,
    )

    apply_btn = widgets.Button(description='Apply Mitigation', button_style='primary')
    approve_btn = widgets.Button(description='Approve All', button_style='success')
    redetect_btn = widgets.Button(description='Re-detect', button_style='info')

    status_label = widgets.Label(value='')
    output_area = widgets.Output()

    def update_display(z_index=None):
        """Update the visualization for current z-plane."""
        if z_index is None:
            z_index = z_slider.value

        with output_area:
            clear_output(wait=True)

            fig, axes = plt.subplots(2, 2, figsize=(12, 10))

            # Original image
            ax1 = axes[0, 0]
            im1 = ax1.imshow(z_stack[z_index], cmap='gray')
            ax1.set_title(f'Original Z-plane {z_index}')
            plt.colorbar(im1, ax=ax1, fraction=0.046)

            # FFT visualization
            ax2 = axes[0, 1]
            detection = report.results.get(z_index, detect_stripe_artifacts(z_stack[z_index]))
            log_power, profile, peaks = get_fft_visualization(z_stack[z_index], detection)

            # Show log power spectrum (center crop for visibility)
            h, w = log_power.shape
            crop = log_power[h//4:3*h//4, w//4:3*w//4]
            ax2.imshow(crop, cmap='viridis', aspect='auto')
            ax2.set_title(f'FFT Power Spectrum\nStripe Score: {detection.stripe_score:.3f} ({detection.severity})')

            # Vertical frequency profile
            ax3 = axes[1, 0]
            ax3.plot(profile, 'b-', linewidth=0.5)
            if peaks:
                peak_values = [profile[min(p, len(profile)-1)] for p in peaks if p < len(profile)]
                valid_peaks = [p for p in peaks if p < len(profile)]
                ax3.plot(valid_peaks, peak_values, 'ro', markersize=8, label='Detected peaks')
            ax3.set_xlabel('Frequency (positive half)')
            ax3.set_ylabel('Normalized Power')
            ax3.set_title('Vertical Frequency Profile')
            ax3.legend()

            # Mitigated image (if available)
            ax4 = axes[1, 1]
            if z_index in current_mitigated:
                im4 = ax4.imshow(current_mitigated[z_index], cmap='gray')
                ax4.set_title(f'Mitigated ({current_decisions.get(z_index, "?")})')
                plt.colorbar(im4, ax=ax4, fraction=0.046)
            else:
                ax4.text(0.5, 0.5, 'No mitigation applied yet',
                        ha='center', va='center', transform=ax4.transAxes)
                ax4.set_title('Mitigated (pending)')

            plt.tight_layout()
            display(fig)
            plt.close(fig)

            # Status
            if detection.has_artifact:
                status = f"⚠️ Artifact detected: {detection.severity} (score={detection.stripe_score:.3f})"
                if detection.frequency_peaks:
                    freqs = ', '.join([f'{f:.3f}' for f in detection.frequency_peaks[:3]])
                    status += f"\nPeak frequencies: {freqs}"
            else:
                status = "✓ No artifact detected in this plane"

            if z_index in current_decisions:
                status += f"\n✓ Decision: {current_decisions[z_index]}"

            status_label.value = status

    def on_z_change(change):
        """Handle z-plane slider change."""
        z_index = change['new']
        detection = report.results.get(z_index, detect_stripe_artifacts(z_stack[z_index]))

        # Update recommended method
        has_neighbors = 0 < z_index < nz - 1
        recommended = get_recommended_method(detection, has_neighbors)
        method_dropdown.value = recommended

        update_display(z_index)

    def on_apply(btn):
        """Apply selected mitigation to current plane."""
        z_index = z_slider.value
        method = method_dropdown.value
        detection = report.results.get(z_index, detect_stripe_artifacts(z_stack[z_index]))

        if method == 'keep_as_is':
            current_decisions[z_index] = 'keep_as_is'
            if z_index in current_mitigated:
                del current_mitigated[z_index]
        else:
            mit_result = mitigate_artifact(
                z_stack[z_index],
                method=method,
                frequencies=detection.frequency_peaks,
                z_stack=z_stack,
                z_index=z_index,
            )
            current_decisions[z_index] = method
            current_mitigated[z_index] = mit_result.image
            print(f"  {mit_result.message}")

        update_display(z_index)

        # Move to next affected plane
        remaining = [p for p in report.affected_planes if p not in current_decisions]
        if remaining:
            z_slider.value = remaining[0]

    def on_approve(btn):
        """Approve all decisions and close."""
        result.decisions = current_decisions.copy()
        result.mitigated_planes = current_mitigated.copy()
        result.skipped_planes = [z for z, m in current_decisions.items() if m == 'keep_as_is']
        result.approved = True

        # Save report if output_dir provided
        if output_dir:
            _save_artifact_report(output_dir, channel_name, report, current_decisions, sensitivity)

        print(f"\n✓ Approved. {len(current_mitigated)} planes mitigated, {len(result.skipped_planes)} kept as-is.")

    def on_redetect(btn):
        """Re-run detection with new sensitivity."""
        nonlocal report
        new_sens = sensitivity_slider.value
        print(f"Re-scanning with sensitivity={new_sens}...")
        report = scan_zstack_for_artifacts(z_stack, sensitivity=new_sens, verbose=False)
        print(f"  {report}")
        update_display()

    # Connect events
    z_slider.observe(on_z_change, names='value')
    apply_btn.on_click(on_apply)
    approve_btn.on_click(on_approve)
    redetect_btn.on_click(on_redetect)

    # Layout
    controls = widgets.VBox([
        widgets.HBox([z_slider, sensitivity_slider, redetect_btn]),
        widgets.HBox([method_dropdown, apply_btn, approve_btn]),
        status_label,
    ])

    # Display
    display(widgets.VBox([controls, output_area]))

    # Initial display
    update_display()

    return result


def _auto_mitigate(
    z_stack: np.ndarray,
    report,
    result: ArtifactTunerResult,
    output_dir: Path | None,
) -> ArtifactTunerResult:
    """Automatic mitigation when widgets unavailable."""
    from ..qc.stripe_mitigation import get_recommended_method, mitigate_artifact

    nz = z_stack.shape[0]

    for z_index in report.affected_planes:
        detection = report.results[z_index]
        has_neighbors = 0 < z_index < nz - 1
        method = get_recommended_method(detection, has_neighbors)

        if method == 'keep_as_is':
            result.decisions[z_index] = 'keep_as_is'
            result.skipped_planes.append(z_index)
        else:
            mit_result = mitigate_artifact(
                z_stack[z_index],
                method=method,
                frequencies=detection.frequency_peaks,
                z_stack=z_stack,
                z_index=z_index,
            )
            result.decisions[z_index] = method
            result.mitigated_planes[z_index] = mit_result.image
            print(f"  Z-plane {z_index}: {mit_result.message}")

    result.approved = True

    if output_dir:
        _save_artifact_report(output_dir, result.channel_name, report, result.decisions, 0.5)

    return result


def _save_artifact_report(
    output_dir: Path,
    channel_name: str,
    report,
    decisions: dict,
    sensitivity: float,
):
    """Save artifact report and metadata to JSON."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    report_data = {
        'scan_timestamp': datetime.now().isoformat(),
        'channel': channel_name,
        'sensitivity': sensitivity,
        'total_planes': report.total_planes,
        'affected_planes': report.affected_planes,
        'worst_severity': report.worst_severity,
        'planes': {}
    }

    for z_index, detection in report.results.items():
        plane_data = {
            'has_artifact': detection.has_artifact,
            'severity': detection.severity,
            'stripe_score': detection.stripe_score,
            'frequency_peaks': detection.frequency_peaks,
        }

        if z_index in decisions:
            plane_data['mitigation'] = decisions[z_index]
            if decisions[z_index] != 'keep_as_is':
                plane_data['mitigated_file'] = f'mitigated/{z_index:03d}.tif'

        report_data['planes'][str(z_index)] = plane_data

    report_path = output_dir / 'artifact_report.json'
    with open(report_path, 'w') as f:
        json.dump(report_data, f, indent=2)

    print(f"  Saved artifact report to {report_path}")


def save_mitigated_stack(
    result: ArtifactTunerResult,
    z_stack: np.ndarray,
    output_dir: Path,
) -> dict:
    """
    Save mitigated z-planes to the mitigated/ subdirectory.

    Parameters
    ----------
    result : ArtifactTunerResult
        Result from artifact_detection_tuner
    z_stack : np.ndarray
        Original z-stack (for reference dtype)
    output_dir : Path
        Base output directory (mitigated/ will be created inside)

    Returns
    -------
    dict
        Mapping of z_index to saved file path
    """
    from skimage.io import imsave

    mitigated_dir = Path(output_dir) / 'mitigated'
    mitigated_dir.mkdir(parents=True, exist_ok=True)

    saved_files = {}

    for z_index, image in result.mitigated_planes.items():
        filename = f'{z_index:03d}.tif'
        filepath = mitigated_dir / filename
        imsave(filepath, image, check_contrast=False)
        saved_files[z_index] = filepath
        print(f"  Saved mitigated z-plane {z_index} to {filepath}")

    return saved_files


def visualize_artifact_comparison(
    original: np.ndarray,
    mitigated: np.ndarray,
    title: str = "",
):
    """
    Display before/after comparison of artifact mitigation.

    Parameters
    ----------
    original : np.ndarray
        Original image with artifact
    mitigated : np.ndarray
        Mitigated image
    title : str
        Plot title
    """
    if not MATPLOTLIB_AVAILABLE:
        print("Matplotlib not available for visualization")
        return

    from ..qc.stripe_artifact import detect_stripe_artifacts

    orig_result = detect_stripe_artifacts(original)
    mit_result = detect_stripe_artifacts(mitigated)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Original
    axes[0].imshow(original, cmap='gray')
    axes[0].set_title(f'Original\nScore: {orig_result.stripe_score:.3f} ({orig_result.severity})')

    # Mitigated
    axes[1].imshow(mitigated, cmap='gray')
    axes[1].set_title(f'Mitigated\nScore: {mit_result.stripe_score:.3f} ({mit_result.severity})')

    # Difference
    diff = original.astype(np.float64) - mitigated.astype(np.float64)
    vmax = np.percentile(np.abs(diff), 99)
    axes[2].imshow(diff, cmap='RdBu', vmin=-vmax, vmax=vmax)
    axes[2].set_title('Difference (removed artifact)')

    if title:
        fig.suptitle(title, fontsize=14)

    plt.tight_layout()
    display(fig)
    plt.close(fig)
