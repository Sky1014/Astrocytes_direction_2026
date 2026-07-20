%% Demo: read the original MATLAB release files
% This script demonstrates the minimal workflow for reading one fish record
% from the NB astroglia direction data release.

script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'astro_functions'));

clear
close all

%% User settings
release_dir = 'F:\Data\Lightsheet\Astrocytes_direction\NB_data_release';
data_name = 'fish_4';
roi_id = 1;
sampling_freq = 6000;

swim_dir = fullfile(release_dir, data_name, 'swimming', [data_name, '_processed']);
image_dir = fullfile(release_dir, data_name, 'imaging', [data_name, '_registered']);

%% Load swim, stimulus, and alignment files
load(fullfile(swim_dir, [data_name, '_right_swim.mat']), 'right_swim');
load(fullfile(swim_dir, [data_name, '_left_swim.mat']), 'left_swim');
load(fullfile(swim_dir, [data_name, '_orient.mat']), 'orient');
load(fullfile(swim_dir, [data_name, '_trial_stage.mat']), 'trial_stage');
load(fullfile(swim_dir, [data_name, '_frames.mat']), 'frames');
load(fullfile(swim_dir, [data_name, '_stages.mat']), 'stages');
load(fullfile(swim_dir, [data_name, '_stages_framed.mat']), 'stages_framed');
load(fullfile(swim_dir, [data_name, '_trials.mat']), 'trials');

fprintf('Loaded %s swim/stimulus vectors: %d samples at %d Hz.\n', ...
    data_name, numel(right_swim), sampling_freq);
fprintf('Loaded %d imaging frames.\n', size(frames, 1));

%% Load ROI metadata and calcium traces
load(fullfile(image_dir, 'cell_info.mat'), 'cell_info');
load(fullfile(image_dir, 'cell_resp_dim_processed.mat'), 'cell_resp_dim');

fprintf('Loading fluorescence matrix from stackf. This can take time for large fish records.\n');
cell_resp_processed = read_LSstack_fast_float_astrodir2026( ...
    fullfile(image_dir, 'cell_resp_processed.stackf'), cell_resp_dim);

fprintf('Fluorescence matrix size: %d ROIs x %d frames.\n', ...
    size(cell_resp_processed, 1), size(cell_resp_processed, 2));

%% Inspect the selected ROI
roi_center_yx = cell_info(roi_id).center;
roi_plane = cell_info(roi_id).slice;
roi_area_px = cell_info(roi_id).area;

fprintf('ROI %d: center [y x] = [%d %d], plane = %d, area = %d pixels.\n', ...
    roi_id, roi_center_yx(1), roi_center_yx(2), roi_plane, roi_area_px);

%% Plot one ROI trace aligned to one visual-motion epoch
% Direction IDs are 0-11 in the manuscript description. MATLAB struct index
% is direction_id + 1.
direction_id = 0;
direction_index = direction_id + 1;
trial_index = 1;

epoch_on_sample = stages(direction_index).onset(trial_index);
epoch_off_sample = stages(direction_index).offset(trial_index);
plot_margin_s = 5;

plot_start_sample = max(1, epoch_on_sample - round(plot_margin_s * sampling_freq));
plot_stop_sample = min(numel(right_swim), epoch_off_sample + round(plot_margin_s * sampling_freq));
sample_window = plot_start_sample:plot_stop_sample;

frame_mask = frames(:, 2) >= plot_start_sample & frames(:, 2) <= plot_stop_sample;
frame_ids = find(frame_mask);
frame_times_s = frames(frame_mask, 2) / sampling_freq;
roi_trace = cell_resp_processed(roi_id, frame_ids);

sample_times_s = sample_window / sampling_freq;
right_swim_window = normalize_trace_for_display(right_swim(sample_window));
left_swim_window = normalize_trace_for_display(left_swim(sample_window));
orient_window = orient(sample_window);
trial_stage_window = trial_stage(sample_window);

figure('Color', 'w', 'Position', [100 100 1100 700]);

subplot(4, 1, 1);
plot(frame_times_s, roi_trace, 'k', 'LineWidth', 1);
ylabel('F/F0 + 1');
title(sprintf('%s ROI %d fluorescence', data_name, roi_id), 'Interpreter', 'none');
box off

subplot(4, 1, 2);
plot(sample_times_s, right_swim_window, 'Color', [0.47 0.76 0.36], 'LineWidth', 0.8);
hold on
plot(sample_times_s, -left_swim_window, 'Color', [0.83 0.30 0.62], 'LineWidth', 0.8);
ylabel('Swim (a.u.)');
legend({'right swim', 'left swim'}, 'Location', 'best');
box off

subplot(4, 1, 3);
plot(sample_times_s, orient_window, 'b', 'LineWidth', 0.8);
ylabel('orient');
box off

subplot(4, 1, 4);
plot(sample_times_s, trial_stage_window, 'r', 'LineWidth', 0.8);
ylabel('trial stage');
xlabel('Time (s)');
box off

for ax = findall(gcf, 'Type', 'axes').'
    xline(ax, epoch_on_sample / sampling_freq, '--', 'Stim on');
    xline(ax, epoch_off_sample / sampling_freq, '--', 'Stim off');
end

%% Local display helper
function y = normalize_trace_for_display(x)
x = double(x(:));
x = x - prctile(x, 70);
x(x < 0) = 0;
scale_value = max(x);
if isempty(scale_value) || ~isfinite(scale_value) || scale_value <= 0
    scale_value = 1;
end
y = x / scale_value;
end
