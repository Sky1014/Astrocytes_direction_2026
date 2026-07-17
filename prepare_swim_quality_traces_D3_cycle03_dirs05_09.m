% Prepare zoomed swim traces for selected direction regions in D3 cycle 3.
%
% This script imports swim-related signals, performs preprocessing, and plots
% two direction-region zoom views from the third orientation cycle.

clear; close all; clc;

%% Input dataset

use_nb_data_release = true;

if use_nb_data_release
    release_dir = 'F:\Data\Lightsheet\Astrocytes_direction\NB_data_release';
    data_name = 'fish_4'; % original data_name: 20240121_2_1
    data_id = string(data_name);
    folder = fullfile(release_dir, data_name, 'swimming', [data_name, '_processed']);
else
    raw_dir = 'F:\Data\Lightsheet\Astrocytes_direction\';
    data_name = '20240121_2_1';
    date = string(data_name(:, 1:8));
    data_id = string(data_name(:, 1:12));

    VR_dir = convertStringsToChars(strcat(raw_dir, date, '\VR\'));
    folder = strcat(VR_dir, data_id, '_processed');
end

SamplingFreq = 6000;

%% Import swim and synchronization signals

load(strcat(folder, '\', data_id, '_ch1.mat'), 'ch1');
load(strcat(folder, '\', data_id, '_ch2.mat'), 'ch2');
load(strcat(folder, '\', data_id, '_orient.mat'), 'orient');
load(strcat(folder, '\', data_id, '_stages.mat'), 'stages');

%% Filter raw swim channels

filtdata1_full = filter_data(ch1);
filtdata2_full = filter_data(ch2);

out_dir = 'D:\WSJ\Mulab\Paper_inbox\astroglia_direction\bilateral_swim';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

swim_sample_window = 1898570:9223790;
swim_sample_window = swim_sample_window(swim_sample_window >= 1 & ...
    swim_sample_window <= min([length(filtdata1_full), length(filtdata2_full), length(orient)]));

swim_baseline_prctile = 70;
swim_scale_prctile = 100;

[swim_ch1_norm, swim_ch1_baseline, swim_ch1_scale] = normalize_swim_trace_01( ...
    filtdata1_full(swim_sample_window), swim_baseline_prctile, swim_scale_prctile);
[swim_ch2_norm, swim_ch2_baseline, swim_ch2_scale] = normalize_swim_trace_01( ...
    filtdata2_full(swim_sample_window), swim_baseline_prctile, swim_scale_prctile);

time_swim_window = swim_sample_window / SamplingFreq;
swim_ch1_mirror_ready = -swim_ch1_norm;
swim_ch2_mirror_ready = swim_ch2_norm;

green_swim = [0.47 0.76 0.36];
pink_swim = [0.83 0.30 0.62];

selected_cycle_number = 3;
selected_direction_orders = [5, 9];
moving_pre_time_s = 5;
subplot_scale_time_s = 5;
out_pdf = fullfile(out_dir, ...
    [char(data_id), '_bilateral_swim_cycle03_direction_regions05_09.pdf']);
out_pdf = Name_File_with_Suffix(out_pdf);

[cycle_start, cycle_end, cycle_trial_id] = find_nth_complete_orientation_cycle( ...
    stages, swim_sample_window(1), swim_sample_window(end), selected_cycle_number);
cycle_direction_events = collect_stage_events(stages, cycle_trial_id, SamplingFreq);

if max(selected_direction_orders) > size(cycle_direction_events, 1)
    error('Cycle %d only has %d direction regions; cannot select region %d.', ...
        selected_cycle_number, size(cycle_direction_events, 1), max(selected_direction_orders));
end

region_specs = struct( ...
    'order_id', {}, ...
    'dir_id', {}, ...
    'moving_onset_sample', {}, ...
    'start_sample', {}, ...
    'end_sample', {}, ...
    'stim_bar_intervals', {}, ...
    'direction_event', {});
for region_idx = 1:numel(selected_direction_orders)
    order_id = selected_direction_orders(region_idx);
    direction_event = cycle_direction_events(order_id, :);
    direction_start = round(direction_event(1) * SamplingFreq);
    region_end = round(direction_event(2) * SamplingFreq);
    moving_onset = find_first_moving_onset_sample( ...
        orient, direction_start, region_end);
    region_start = max(swim_sample_window(1), ...
        moving_onset - round(moving_pre_time_s * SamplingFreq));

    region_specs(region_idx).order_id = order_id;
    region_specs(region_idx).dir_id = direction_event(3);
    region_specs(region_idx).moving_onset_sample = moving_onset;
    region_specs(region_idx).start_sample = region_start;
    region_specs(region_idx).end_sample = region_end;
    region_specs(region_idx).stim_bar_intervals = collect_stim_on_intervals( ...
        orient, stages, region_start, region_end, SamplingFreq);
    region_specs(region_idx).direction_event = direction_event;
end

plot_bilateral_swim_direction_regions(time_swim_window, ...
    swim_ch1_mirror_ready, swim_ch2_mirror_ready, swim_sample_window, ...
    green_swim, pink_swim, out_pdf, region_specs, ...
    make_direction_colors(), subplot_scale_time_s);

fprintf('Prepared and plotted swim quality traces for %s.\n', data_id);
fprintf('Full window samples: %d to %d (%.3f to %.3f s).\n', ...
    swim_sample_window(1), swim_sample_window(end), time_swim_window(1), time_swim_window(end));
fprintf('Zoom cycle: trial %d, samples %d to %d (%.3f to %.3f s).\n', ...
    cycle_trial_id, cycle_start, cycle_end, cycle_start / SamplingFreq, cycle_end / SamplingFreq);
for region_idx = 1:numel(region_specs)
    fprintf('Direction region order %d: dir %d, moving onset %d (%.3f s), plot samples %d to %d (%.3f to %.3f s).\n', ...
        region_specs(region_idx).order_id, region_specs(region_idx).dir_id, ...
        region_specs(region_idx).moving_onset_sample, ...
        region_specs(region_idx).moving_onset_sample / SamplingFreq, ...
        region_specs(region_idx).start_sample, region_specs(region_idx).end_sample, ...
        region_specs(region_idx).start_sample / SamplingFreq, ...
        region_specs(region_idx).end_sample / SamplingFreq);
end
fprintf('CH1 baseline/scale: %.6g / %.6g\n', swim_ch1_baseline, swim_ch1_scale);
fprintf('CH2 baseline/scale: %.6g / %.6g\n', swim_ch2_baseline, swim_ch2_scale);


function plot_bilateral_swim_direction_regions(time_sec, ch1_up, ch2_down, ...
    sample_window, green_swim, pink_swim, out_pdf, region_specs, ...
    direction_colors, scale_time_s)
fig = figure('Color', 'w', ...
    'Units', 'pixels', ...
    'Position', [250 250 1180 360], ...
    'Visible', 'off');
n_regions = numel(region_specs);
axes_left_first = 0.07;
axes_bottom = 0.22;
axes_width = 0.40;
axes_height = 0.62;
axes_gap = 0.08;

for region_idx = 1:n_regions
    axes_left = axes_left_first + (region_idx - 1) * (axes_width + axes_gap);
    ax = axes(fig, 'Position', [axes_left axes_bottom axes_width axes_height]);
    region_mask = sample_window >= region_specs(region_idx).start_sample & ...
        sample_window <= region_specs(region_idx).end_sample;
    plot_bilateral_swim_axis(ax, time_sec(region_mask), ...
        ch1_up(region_mask), ch2_down(region_mask), green_swim, pink_swim, ...
        region_specs(region_idx).stim_bar_intervals, ...
        region_specs(region_idx).direction_event, direction_colors, scale_time_s);
end

exportgraphics(fig, out_pdf, ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');
close(fig);
end


function plot_bilateral_swim_axis(ax, time_sec, ch1_up, ch2_down, green_swim, pink_swim, ...
    stim_bar_intervals, direction_events, direction_colors, scale_time_s)
hold(ax, 'on');

plot_stride = max(1, floor(numel(time_sec) / 80000));
plot_idx = 1:plot_stride:numel(time_sec);

has_stim_bar_intervals = ~isempty(stim_bar_intervals);
if has_stim_bar_intervals
    y_patch_bottom = -1.05;
    y_patch_top = 1.05;
    for k = 1:size(stim_bar_intervals, 1)
        patch(ax, ...
            [stim_bar_intervals(k, 1), stim_bar_intervals(k, 2), ...
             stim_bar_intervals(k, 2), stim_bar_intervals(k, 1)], ...
            [y_patch_bottom, y_patch_bottom, y_patch_top, y_patch_top], ...
            [0.2118 0.6118 0.8471], ...
            'FaceAlpha', 0.07, ...
            'EdgeColor', 'none');
    end
end

plot(ax, time_sec(plot_idx), ch1_up(plot_idx), ...
    'Color', green_swim, ...
    'LineWidth', 0.85);
plot(ax, time_sec(plot_idx), ch2_down(plot_idx), ...
    'Color', pink_swim, ...
    'LineWidth', 0.85);
plot(ax, [time_sec(1), time_sec(end)], [0, 0], '-', ...
    'Color', [0.55 0.55 0.55], ...
    'LineWidth', 0.6);

has_direction_events = ~isempty(direction_events);
if has_direction_events
    arrow_y = -1.20;
    arrow_dx = max(2.0, 0.020 * (time_sec(end) - time_sec(1)));
    arrow_dy = 0.12;
    for k = 1:size(direction_events, 1)
        dir_id = direction_events(k, 3);
        mid_s = mean(direction_events(k, 1:2));
        theta = -pi/2 - (dir_id - 1) * pi / 6;
        quiver(ax, mid_s, arrow_y, cos(theta) * arrow_dx, sin(theta) * arrow_dy, ...
            0, ...
            'Color', direction_colors(dir_id, :), ...
            'LineWidth', 1.4, ...
            'MaxHeadSize', 1.2);
    end
    ylim(ax, [-1.35 1.08]);
else
    ylim(ax, [-1.08 1.08]);
end

xlim(ax, [time_sec(1), time_sec(end)]);
draw_time_scale_bar(ax, scale_time_s);
style_axes(ax);
end


function [cycle_start, cycle_end, trial_id] = find_nth_complete_orientation_cycle( ...
    stages, window_start, window_end, nth_cycle)
n_trials = min(arrayfun(@(s) numel(s.onset), stages));
overlap_trials = [];

for t = 1:n_trials
    cycle_start_candidate = min(arrayfun(@(s) s.onset(t), stages));
    cycle_end_candidate = max(arrayfun(@(s) s.offset(t), stages));
    if cycle_end_candidate >= window_start && cycle_start_candidate <= window_end
        overlap_trials(end + 1) = t; %#ok<AGROW>
    end
end

if numel(overlap_trials) < nth_cycle
    error('Only %d orientation cycles overlap with the swim window; cannot select cycle %d.', ...
        numel(overlap_trials), nth_cycle);
end

trial_id = overlap_trials(nth_cycle);
cycle_start = min(arrayfun(@(s) s.onset(trial_id), stages));
cycle_end = max(arrayfun(@(s) s.offset(trial_id), stages));
end


function stage_events = collect_stage_events(stages, trial_id, sampling_freq)
stage_events = zeros(numel(stages), 3);
for dir_id = 1:numel(stages)
    stage_events(dir_id, 1) = stages(dir_id).onset(trial_id) / sampling_freq;
    stage_events(dir_id, 2) = stages(dir_id).offset(trial_id) / sampling_freq;
    stage_events(dir_id, 3) = dir_id;
end
stage_events = sortrows(stage_events, 1);
end


function moving_onset = find_first_moving_onset_sample(orient, window_start, window_end)
moving = abs(orient) > 0;
moving_edges = diff([false, moving(:).', false]);
moving_onsets = find(moving_edges == 1);
moving_offsets = find(moving_edges == -1) - 1;

overlap_mask = moving_offsets >= window_start & moving_onsets <= window_end;
if any(overlap_mask)
    overlapping_onsets = moving_onsets(overlap_mask);
    moving_onset = max(overlapping_onsets(1), window_start);
else
    moving_onset = window_start;
end
end


function stim_bar_intervals = collect_stim_on_intervals( ...
    orient, stages, window_start, window_end, sampling_freq)
stim_bar_intervals = [];

moving = abs(orient) > 0;
moving_edges = diff([false, moving(:).', false]);
moving_onsets = find(moving_edges == 1);
moving_offsets = find(moving_edges == -1) - 1;

for moving_id = 1:length(moving_onsets)
    onset_sample = moving_onsets(moving_id);
    offset_sample = moving_offsets(moving_id);
    if offset_sample >= window_start && onset_sample <= window_end
        clipped_onset = max(onset_sample, window_start);
        clipped_offset = min(offset_sample, window_end);
        stim_bar_intervals(end + 1, :) = [ ...
            clipped_onset / sampling_freq, ...
            clipped_offset / sampling_freq]; %#ok<AGROW>
    end
end

moving_bar_duration = [];
if ~isempty(moving_onsets)
    moving_bar_duration = median((moving_offsets - moving_onsets + 1) / sampling_freq);
end
if isempty(moving_bar_duration) || moving_bar_duration <= 0
    moving_bar_duration = 10;
end

initial_dir_id = 1;
if initial_dir_id <= length(stages)
    for trial_id = 1:length(stages(initial_dir_id).onset)
        onset_s = stages(initial_dir_id).onset(trial_id) / sampling_freq;
        offset_s = onset_s + moving_bar_duration;
        if offset_s >= window_start / sampling_freq && onset_s <= window_end / sampling_freq
            clipped_onset_s = max(onset_s, window_start / sampling_freq);
            clipped_offset_s = min(offset_s, window_end / sampling_freq);
            stim_bar_intervals(end + 1, :) = [clipped_onset_s, clipped_offset_s]; %#ok<AGROW>
        end
    end
end

if ~isempty(stim_bar_intervals)
    stim_bar_intervals = sortrows(stim_bar_intervals, 1);
    stim_bar_intervals = unique(round(stim_bar_intervals * sampling_freq) / sampling_freq, 'rows');
end
end


function colors = make_direction_colors()
colors = hsv(12);
colors_hsv = rgb2hsv(colors);
colors_hsv(:, 2) = colors_hsv(:, 2) * 0.65;
colors_hsv(:, 3) = colors_hsv(:, 3) * 0.85;
colors = hsv2rgb(colors_hsv);
end


function style_axes(ax)
set(ax, ...
    'Visible', 'off', ...
    'Box', 'off', ...
    'XTick', [], ...
    'YTick', [], ...
    'FontName', 'Arial', ...
    'FontSize', 12);
grid(ax, 'off');
end


function draw_time_scale_bar(ax, scale_time_s)
x_limits = xlim(ax);
y_limits = ylim(ax);
x_span = x_limits(2) - x_limits(1);
y_span = y_limits(2) - y_limits(1);

x2 = x_limits(2) - 0.04 * x_span;
x1 = x2 - scale_time_s;
y = y_limits(2) - 0.08 * y_span;
text_gap = 0.035 * y_span;

plot(ax, [x1, x2], [y, y], 'k-', ...
    'LineWidth', 1.8, ...
    'Clipping', 'off');
text(ax, (x1 + x2) / 2, y + text_gap, sprintf('%g s', scale_time_s), ...
    'Color', 'k', ...
    'FontName', 'Arial', ...
    'FontSize', 13, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'Clipping', 'off');
end


function [y, baseline, scale_value] = normalize_swim_trace_01(x, baseline_prctile, scale_prctile)
x = double(x(:));
finite_mask = isfinite(x);
if ~any(finite_mask)
    y = zeros(size(x));
    baseline = 0;
    scale_value = 1;
    return;
end

baseline = local_percentile(x(finite_mask), baseline_prctile);
y = x - baseline;
y(~isfinite(y)) = 0;
y(y < 0) = 0;

scale_value = local_percentile(y(y > 0), scale_prctile);
if isempty(scale_value) || ~isfinite(scale_value) || scale_value <= 0
    scale_value = max(y);
end
if isempty(scale_value) || ~isfinite(scale_value) || scale_value <= 0
    scale_value = 1;
end

y = min(y / scale_value, 1);
end


function p = local_percentile(x, percentile_value)
x = sort(x(:));
x = x(isfinite(x));
if isempty(x)
    p = [];
    return;
end

idx = round(percentile_value / 100 * numel(x));
idx = max(1, min(numel(x), idx));
p = x(idx);
end
