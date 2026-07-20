"""Demo: read the lightweight NWB supplement for the NB astroglia dataset.

This script demonstrates the minimal workflow for reading one fish-level NWB
file, inspecting ROI metadata, and plotting one ROI trace aligned with swim and
stimulus variables.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from pynwb import NWBHDF5IO


def normalize_trace_for_display(x: np.ndarray) -> np.ndarray:
    """Baseline-subtract and scale one trace for compact plotting."""
    x = np.asarray(x, dtype=float).reshape(-1)
    x = x - np.nanpercentile(x, 70)
    x[x < 0] = 0
    scale_value = np.nanmax(x) if x.size else 1
    if not np.isfinite(scale_value) or scale_value <= 0:
        scale_value = 1
    return x / scale_value


def main() -> None:
    release_dir = Path(r"F:\Data\Lightsheet\Astrocytes_direction\NB_data_release")
    fish_id = "fish_4"
    roi_id = 0  # Python/NWB uses 0-based ROI indexing.
    sampling_freq = 6000.0

    nwb_path = release_dir / "nwb" / f"{fish_id}.nwb"
    if not nwb_path.exists():
        raise FileNotFoundError(f"NWB file not found: {nwb_path}")

    with NWBHDF5IO(str(nwb_path), "r", load_namespaces=True) as io:
        nwb = io.read()

        fluorescence = nwb.processing["ophys"]["baseline_offset_normalized_fluorescence"]
        right_swim = nwb.acquisition["right_swim"]
        left_swim = nwb.acquisition["left_swim"]
        orient = nwb.stimulus["orient"]
        trial_stage = nwb.stimulus["trial_stage"]
        roi_metadata = nwb.scratch["astrocytic_roi_metadata"].data
        roi_metadata_columns = nwb.scratch["astrocytic_roi_metadata_columns"].data
        epochs = nwb.intervals["visual_motion_stimulus_epochs"].to_dataframe()

        print(f"Session ID: {nwb.session_id}")
        print(f"Subject: {nwb.subject.subject_id}, {nwb.subject.species}, {nwb.subject.age}")
        print(f"Fluorescence shape (frame x ROI): {fluorescence.data.shape}")
        print(f"Right swim samples: {right_swim.data.shape[0]}")
        print(f"Left swim samples: {left_swim.data.shape[0]}")
        print(f"Stimulus epochs: {len(epochs)}")
        print(f"ROI metadata columns: {roi_metadata_columns}")
        print(f"Selected ROI metadata: {roi_metadata[roi_id, :]}")

        # Use the first visual-motion epoch as a compact alignment example.
        epoch = epochs.iloc[0]
        start_sample = int(epoch["ephys_onset_sample_1based"])
        stop_sample = int(epoch["ephys_offset_sample_1based"])
        margin_samples = int(5 * sampling_freq)
        plot_start = max(1, start_sample - margin_samples)
        plot_stop = min(right_swim.data.shape[0], stop_sample + margin_samples)

        # Convert 1-based MATLAB/ephys sample indices to 0-based Python slices.
        sample_slice = slice(plot_start - 1, plot_stop)
        sample_times_s = np.arange(plot_start, plot_stop + 1) / sampling_freq

        frame_timestamps = fluorescence.timestamps[:]
        frame_mask = (frame_timestamps >= plot_start / sampling_freq) & (
            frame_timestamps <= plot_stop / sampling_freq
        )
        frame_indices = np.flatnonzero(frame_mask)
        frame_times_s = frame_timestamps[frame_indices]
        roi_trace = fluorescence.data[frame_indices, roi_id]

        right_swim_window = normalize_trace_for_display(right_swim.data[sample_slice])
        left_swim_window = normalize_trace_for_display(left_swim.data[sample_slice])
        orient_window = orient.data[sample_slice]
        trial_stage_window = trial_stage.data[sample_slice]

    fig, axes = plt.subplots(4, 1, figsize=(11, 7), sharex=False)

    axes[0].plot(frame_times_s, roi_trace, color="black", linewidth=1)
    axes[0].set_ylabel("F/F0 + 1")
    axes[0].set_title(f"{fish_id} ROI {roi_id} fluorescence")

    axes[1].plot(sample_times_s, right_swim_window, color=(0.47, 0.76, 0.36), linewidth=0.8)
    axes[1].plot(sample_times_s, -left_swim_window, color=(0.83, 0.30, 0.62), linewidth=0.8)
    axes[1].set_ylabel("Swim (a.u.)")
    axes[1].legend(["right swim", "left swim"], loc="best")

    axes[2].plot(sample_times_s, orient_window, color="blue", linewidth=0.8)
    axes[2].set_ylabel("orient")

    axes[3].plot(sample_times_s, trial_stage_window, color="red", linewidth=0.8)
    axes[3].set_ylabel("trial stage")
    axes[3].set_xlabel("Time (s)")

    for ax in axes:
        ax.axvline(start_sample / sampling_freq, color="gray", linestyle="--", linewidth=0.8)
        ax.axvline(stop_sample / sampling_freq, color="gray", linestyle="--", linewidth=0.8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
