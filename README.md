# Code Description

This folder contains MATLAB scripts used to generate example figures and quality-control plots for the astrocyte direction dataset.

## `D3_NB_Figure_withSwim.m`

Generates an example figure combining directional visual stimulation, astrocytic calcium activity, and swim signals. The script loads `right_swim/left_swim`, `orient`, `frames`, `stages`, `cell_resp_processed.stackf`, and `cell_info.mat` from one fish, then aligns swim signals, visual stimuli, and calcium imaging frames to a common time axis. It identifies direction-related astrocytic ROIs using 12 direction-specific stimulus kernels, plots direction-colored average calcium traces, and overlays stimulus bars, direction markers, and swim traces. It also exports a separate figure showing the 12 direction kernels/regressors.

## `NB_direction_cell_fluorescence_map_20260704.m`

Plots the spatial distribution of direction-preferring astrocytic ROIs. The script loads the calcium activity matrix, ROI locations, stimulus timing, and direction information, then estimates each ROI's preferred direction using either a regressor-based or trial-peak-based method. Direction-preferring ROIs are overlaid on an average brain image or registered template, with colors indicating preferred visual-motion direction.

## `prepare_swim_quality_traces_D3.m`

Generates bilateral swim recording-quality overview plots. The script loads `right_swim/left_swim`, `orient`, and `stages`, selects a predefined long recording window, and applies display-oriented percentile baseline subtraction and scaling to the swim traces. It exports one PDF showing a long-window bilateral swim overview and another PDF showing a zoomed view of the third orientation cycle, with stimulus shading and direction markers.

## `prepare_swim_quality_traces_D3_cycle03_dirs05_09.m`

Generates more focused zoom-in plots for two selected direction regions within the third orientation cycle. The current script selects the 5th and 9th direction regions. It loads and normalizes the swim traces, identifies the third complete orientation cycle, and plots local bilateral swim traces around the selected direction stimuli. The output is a PDF with two side-by-side subplots for inspecting swim recording quality near those stimulus periods.

## `plot_motion_param_depth_summary.m`

Creates depth-resolved summary plots of motion drift from `motion_param.mat`. The script reads segmented-grid motion estimates for each imaging plane, calculates mean and SEM values for XY displacement and Z displacement, and exports two PDF summaries: one for XY motion and one for Z motion.

## `plot_plane12_xy_z_tcourse_shaded_seconds.m`

Plots XY and Z motion time courses for imaging plane 12. The script reads `move_tcourse_for_plot.mat` and uses `frames.mat` to convert frame indices to seconds. It exports shaded time-course PDFs for Z displacement and XY displacement over recording time, with the shaded region representing either SD or SEM. This script is used to show how drift changes over time in one selected imaging plane.

## `plot_plane12_motion_patches.m`

Visualizes segmented-grid motion patches on the average image of imaging plane 12. The script reads `ave.tif` and `motion_param.mat`, then generates two motion quality-control figures: one showing Z displacement as colored square patches and another showing XY displacement as red arrows. It exports two PDFs that show the spatial distribution of estimated motion drift on the imaging plane.
