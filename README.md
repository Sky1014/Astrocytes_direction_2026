# Code Description

This repository provides MATLAB example code for reading, aligning, quality-checking, and analyzing processed whole-brain astrocytic calcium imaging data, frame-aligned visual-motion stimulus traces, bilateral fictive-swimming recordings, and ROI metadata associated with direction-dependent astroglial responses in larval zebrafish.

## Repository structure

- Root directory: main MATLAB scripts for figure generation, swim quality-control plots, and motion quality-control plots.
- `astro_functions/`: bundled helper functions used by the main scripts. These helper functions use the `_astrodir2026` suffix to avoid name collisions with similarly named functions in other MATLAB projects.
- `README.md`: this file.

Each main script adds `astro_functions/` to the MATLAB path at runtime, so the scripts can be run from the repository root without adding external project-specific helper folders.

## Main scripts

### `D3_NB_Figure_withSwim.m`

Generates an example figure combining directional visual stimulation, astrocytic calcium activity, and swim signals. The script loads swim traces, visual stimulus timing, imaging frame timing, ROI calcium traces, and ROI metadata from one fish record. It identifies direction-related astrocytic ROIs using 12 direction-specific stimulus kernels, plots direction-colored average calcium traces, and overlays stimulus bars, direction markers, and swim traces.

### `NB_direction_cell_fluorescence_map_20260704.m`

Plots the spatial distribution of direction-preferring astrocytic ROIs. The script estimates each ROI's preferred direction using either a regressor-based or trial-peak-based method and overlays selected ROIs on an average brain image or registered template.

### `prepare_swim_quality_traces_D3.m`

Generates bilateral swim recording-quality overview plots. The script loads swim traces, visual stimulus timing, and stage information, applies display-oriented percentile baseline subtraction and scaling, and exports overview and zoomed swim quality-control figures.

### `prepare_swim_quality_traces_D3_cycle03_dirs05_09.m`

Generates focused zoom-in plots for two selected direction regions within the third orientation cycle. The output is a two-panel PDF for inspecting bilateral swim recording quality around selected stimulus periods.

### `plot_motion_param_depth_summary.m`

Creates depth-resolved summary plots of motion drift from `motion_param.mat`. The script calculates mean and SEM values for XY displacement and Z displacement across segmented-grid regions and exports PDF summaries.

### `plot_plane12_xy_z_tcourse_shaded_seconds.m`

Plots XY and Z motion time courses for imaging plane 12. The script reads saved motion time-course data and uses `frames.mat` to convert frame indices to seconds.

### `plot_plane12_motion_patches.m`

Visualizes segmented-grid motion patches on the average image of imaging plane 12. The script generates quality-control figures showing Z displacement as colored patches and XY displacement as arrows.

## Helper functions

The `astro_functions/` folder contains local copies of the custom functions required by the scripts, including file naming, swim preprocessing, TIFF/NRRD reading, light-sheet stack reading, direction-regressor construction, and colormap generation. The bundled `read_LSstack_fast_float_mex64_astrodir2026.mexw64` file is required by `read_LSstack_fast_float_astrodir2026.m`, and `slanCM_Data.mat` is required by `slanCM_astrodir2026.m`.
