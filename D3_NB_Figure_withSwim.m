%% import ephys data
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'astro_functions'));

disp('Clearing workspace.');
clear all
close all

colormap_glia = [];
colormap_glia = colormap(slanCM_astrodir2026('gist_rainbow'));
close all

tic
use_nb_data_release = true;

x_window = [542,1040];

if use_nb_data_release
    release_dir = 'F:\Data\Lightsheet\Astrocytes_direction\NB_data_release';
    data_name = 'fish_4';
    data_id = string(data_name);
    folder = fullfile(release_dir, data_name, 'swimming', [data_name, '_processed']);
    image_dir = fullfile(release_dir, data_name, 'imaging', [data_name, '_registered']);
else
    raw_dir = 'F:\Data\Lightsheet\Astrocytes_direction\';
    data_name = '20240121_2_1';
    date = string(data_name(:,1:8));
    data_id = string(data_name(:,1:12));

    VR_dir = convertStringsToChars(strcat(raw_dir,date,'\VR\'));
    Image_dir = convertStringsToChars(strcat(raw_dir,date,'\Volume\'));
    folder = strcat(VR_dir,data_id,'_processed');
    image_dir = convertStringsToChars(strcat(Image_dir,string(data_id),'_registered'));
end

if use_nb_data_release
    load(strcat(folder,'\',data_id,'_right_swim.mat'), 'right_swim');
    load(strcat(folder,'\',data_id,'_left_swim.mat'), 'left_swim');
    filtdata1 = right_swim;
    filtdata2 = left_swim;
else
    load(strcat(folder,'\',data_id,'_ch1.mat'));
    load(strcat(folder,'\',data_id,'_ch2.mat'));
    filtdata1 = filter_data_astrodir2026(ch1);
    filtdata2 = filter_data_astrodir2026(ch2);
end

load(strcat(folder,'\',data_id,'_orient.mat'));
load(strcat(folder,'\',data_id,'_frames.mat'));
load(strcat(folder,'\',data_id,'_trials.mat'));
load(strcat(folder,'\',data_id,'_stages.mat'));
load(strcat(folder,'\',data_id,'_stages_framed.mat'));

disp("Loading imaging data...")
cell_resp_dim_processed = [];
load(fullfile(image_dir,'\cell_resp_dim_processed.mat'));
cell_resp_processed = [];
cell_resp_processed = read_LSstack_fast_float_astrodir2026(fullfile(image_dir,'\cell_resp_processed.stackf'),cell_resp_dim);
disp('cell_resp_processed loaded');
cell_info = [];
if exist(fullfile(image_dir,'\cell_loc.csv'))
    load(fullfile(image_dir,'\cell_info_processed.mat'));
    disp('Registered cell locations loaded.');
    cells = [];
    for k = 1: cell_resp_dim(1)
        cells(k,1) = cell_info(k).center_aligned(2);
        cells(k,2) = cell_info(k).center_aligned(1);
        cells(k,3) = cell_info(k).center_aligned(3);
    end
else
    load(fullfile(image_dir,'\cell_info.mat'));
    disp('Unregistered cell locations loaded.');
end

m = 1;
n = cell_resp_dim(2);

SamplingFreq = 6000;
frame_start = max(1, frames(1,2));
frame_end = min([frames(end,2), length(filtdata1), length(filtdata2), length(orient)]);
swim_samples = frame_start:frame_end;
filtdata1 = filtdata1(swim_samples);
filtdata2 = filtdata2(swim_samples);
orient = orient(swim_samples);
x_scale_swim = swim_samples / SamplingFreq;
x_scale_frames = (frames(:,2)/SamplingFreq).';

swim_baseline_prctile = 70;
swim_scale_prctile = 100;
swim_ch1_norm_full = normalize_swim_trace_01(filtdata1, ...
    swim_baseline_prctile, swim_scale_prctile);
swim_ch2_norm_full = normalize_swim_trace_01(filtdata2, ...
    swim_baseline_prctile, swim_scale_prctile);
swim_ch1_mirror_ready_full = -swim_ch1_norm_full;
swim_ch2_mirror_ready_full = swim_ch2_norm_full;

num_directions = 12;
trace_spacing = 0.14;

row_ids = round(linspace(1, size(colormap_glia,1), num_directions));
row_ids(row_ids < 1) = 1;
row_ids(row_ids > size(colormap_glia,1)) = size(colormap_glia,1);
trace_colors = colormap_glia(row_ids,:);
trace_hsv = rgb2hsv(trace_colors);
trace_hsv(:,2) = trace_hsv(:,2) * 0.60;
trace_hsv(:,3) = trace_hsv(:,3) * 0.82;
trace_colors = hsv2rgb(trace_hsv);

%%
shuffled_trials = 0;
Exp = [];
Exp(1,1) = 1;
if shuffled_trials == 0
    Exp(2,1) = stages(1).onset(1);
    Exp(3,1) = stages(11).offset(end);
end

rho_standard = 0.1;
mode1_trace_cache_dir = 'D:\WSJ\Mulab\Paper_inbox\astroglia_direction\colorful_calcium_traces';
if ~exist(mode1_trace_cache_dir, 'dir')
    mkdir(mode1_trace_cache_dir);
end
mode1_trace_cache_path = fullfile(mode1_trace_cache_dir, ...
    [char(data_name), '_mode1_trace_cache.mat']);
mode1_trace_cache_version = 1;
use_mode1_trace_cache = false;

if exist(mode1_trace_cache_path, 'file')
    mode1_trace_cache = load(mode1_trace_cache_path);
    cache_has_required_fields = isfield(mode1_trace_cache, 'mode1_trace_mean_all') && ...
        isfield(mode1_trace_cache, 'mode1_cell_index_all') && ...
        isfield(mode1_trace_cache, 'cache_data_name') && ...
        isfield(mode1_trace_cache, 'cache_version') && ...
        isfield(mode1_trace_cache, 'cache_exp') && ...
        isfield(mode1_trace_cache, 'cache_rho_standard') && ...
        isfield(mode1_trace_cache, 'cache_num_directions');
    if cache_has_required_fields && ...
            mode1_trace_cache.cache_version == mode1_trace_cache_version && ...
            strcmp(char(mode1_trace_cache.cache_data_name), char(data_name)) && ...
            mode1_trace_cache.cache_num_directions == num_directions && ...
            isequal(mode1_trace_cache.cache_exp, Exp) && ...
            isequal(mode1_trace_cache.cache_rho_standard, rho_standard) && ...
            size(mode1_trace_cache.mode1_trace_mean_all, 1) == num_directions && ...
            size(mode1_trace_cache.mode1_trace_mean_all, 2) == numel(x_scale_frames)
        mode1_trace_mean_all = mode1_trace_cache.mode1_trace_mean_all;
        mode1_cell_index_all = mode1_trace_cache.mode1_cell_index_all;
        use_mode1_trace_cache = true;
        fprintf('Loaded cached Mode1 traces from %s.\n', mode1_trace_cache_path);
    else
        fprintf('Ignoring stale Mode1 trace cache: %s.\n', mode1_trace_cache_path);
    end
end

if ~use_mode1_trace_cache
    mode1_trace_mean_all = nan(num_directions, numel(x_scale_frames));
    mode1_cell_index_all = cell(num_directions, 1);
end

%% get mode1 timepoints for example

main_fig = figure('Color','w','Position',[100 100 980 540]);
hold on
set(gca,'Color','w','XColor','k','YColor','k','Box','off','TickDir','out');

Mode1_ = [];

for m = 1:12
    if use_mode1_trace_cache
        temp = mode1_trace_mean_all(m, :);
    else
        Mode1_ = [];
        Mode1_(2,:) = stages(m).onset;
        Mode1_(3,:) = stages(m).offset;
        Mode1_(1,:) = 1:size(Mode1_,2);

        frame_qualified = [];
        frame_qualified = frames(find(frames(1:end,2) >= Exp(2,1) & frames(1:end,2) <= Exp(3,end)),2);
        cell_resp_qualified = [];
        cell_resp_qualified = cell_resp_processed(:,find(frames(1:end,2) >= Exp(2,1) & frames(1:end,2) <= Exp(3,end)));
        frame_scale_qualified = (frame_qualified - frame_qualified(1,1) + 1)/6000;

        Mode1_square_kernel = [];
        Mode1_square_kernel = make_full_square_kernel_astrodir2026(Mode1_, Exp);
        Mode1_cells = [];
        Mode1_cells = make_customized_regressor_astrodir2026(Mode1_square_kernel, frame_scale_qualified, frame_qualified, cell_resp_qualified, rho_standard);

        temp = [];
        temp = mean(cell_resp_processed(Mode1_cells.cell_index,:))-1;
        mode1_trace_mean_all(m, :) = temp;
        mode1_cell_index_all{m} = Mode1_cells.cell_index;
    end
    plot(x_scale_frames, temp - m * trace_spacing, 'DisplayName', num2str(m), ...
        'Color', trace_colors(m,:), 'LineWidth', 1.6);
    xlim(x_window);
end

if ~use_mode1_trace_cache
    cache_data_name = data_name;
    cache_version = mode1_trace_cache_version;
    cache_exp = Exp;
    cache_rho_standard = rho_standard;
    cache_num_directions = num_directions;
    mode1_trace_save_path = Name_File_with_Suffix_astrodir2026(mode1_trace_cache_path);
    save(mode1_trace_save_path, 'mode1_trace_mean_all', ...
        'mode1_cell_index_all', 'cache_data_name', 'cache_version', ...
        'cache_exp', 'cache_rho_standard', 'cache_num_directions');
    fprintf('Saved Mode1 trace cache to %s.\n', mode1_trace_save_path);
end

%% add orient-based moving stimulus bars, swim traces, and direction lines

stim_events = [];
for dir_id = 1:num_directions
    if dir_id <= length(stages)
        for trial_id = 1:length(stages(dir_id).onset)
            onset_s = stages(dir_id).onset(trial_id) / SamplingFreq;
            offset_s = stages(dir_id).offset(trial_id) / SamplingFreq;
            if offset_s >= x_window(1) && onset_s <= x_window(2)
                stim_events(end+1,:) = [onset_s, offset_s, dir_id];
            end
        end
    end
end

if ~isempty(stim_events)
    stim_events = sortrows(stim_events, 1);
end

moving_bar_bottom = -trace_spacing * (num_directions + 1.80);
moving_bar_top = 0.20;
swim_base = -trace_spacing * (num_directions + 3.15);
swim_height = 0.26;
direction_line_y = swim_base - 0.18;
direction_line_dx = 9;
direction_line_dy = 0.070;
ds = 50;

stim_bar_intervals = [];
moving = abs(orient) > 0;
moving_edges = diff([false, moving(:).', false]);
moving_onsets = find(moving_edges == 1);
moving_offsets = find(moving_edges == -1) - 1;

for moving_id = 1:length(moving_onsets)
    onset_s = x_scale_swim(moving_onsets(moving_id));
    offset_s = x_scale_swim(moving_offsets(moving_id));
    if offset_s >= x_window(1) && onset_s <= x_window(2)
        stim_bar_intervals(end+1,:) = [onset_s, offset_s];
    end
end

moving_bar_duration = [];
if ~isempty(moving_onsets)
    moving_bar_duration = median((moving_offsets - moving_onsets + 1) / SamplingFreq);
end
if isempty(moving_bar_duration) || moving_bar_duration <= 0
    moving_bar_duration = 10;
end

initial_dir_id = 1;
if initial_dir_id <= length(stages)
    for trial_id = 1:length(stages(initial_dir_id).onset)
        onset_s = stages(initial_dir_id).onset(trial_id) / SamplingFreq;
        offset_s = onset_s + moving_bar_duration;
        if offset_s >= x_window(1) && onset_s <= x_window(2)
            stim_bar_intervals(end+1,:) = [onset_s, offset_s];
        end
    end
end

stim_bar_handles = [];
if ~isempty(stim_bar_intervals)
    stim_bar_intervals = sortrows(stim_bar_intervals, 1);
    stim_bar_intervals = unique(round(stim_bar_intervals * SamplingFreq) / SamplingFreq, 'rows');

    for interval_id = 1:size(stim_bar_intervals,1)
        onset_s = stim_bar_intervals(interval_id,1);
        offset_s = stim_bar_intervals(interval_id,2);
        stim_bar_handles(end+1) = patch([onset_s offset_s offset_s onset_s], ...
            [moving_bar_bottom moving_bar_bottom moving_bar_top moving_bar_top], ...
            [0.2118 0.6118 0.8471], ...
            'FaceAlpha', 0.07, ...
            'EdgeColor', 'none');
    end
    uistack(stim_bar_handles, 'bottom');
end

swim_mask = x_scale_swim >= x_window(1) & x_scale_swim <= x_window(2);
swim_idx = find(swim_mask);
if ~isempty(swim_idx)
    swim_ch1 = swim_ch1_mirror_ready_full(swim_idx);
    swim_ch2 = swim_ch2_mirror_ready_full(swim_idx);
    swim_display_scale = max(abs([swim_ch1(:); swim_ch2(:)]));
    if isempty(swim_display_scale) || ~isfinite(swim_display_scale) || swim_display_scale <= 0
        swim_display_scale = 1;
    end
    swim_ch1 = swim_ch1 / swim_display_scale * swim_height;
    swim_ch2 = swim_ch2 / swim_display_scale * swim_height;

    x_swim_plot = x_scale_swim(swim_idx);
    swim_plot_stride = max(1, floor(numel(x_swim_plot) / 80000));
    swim_plot_idx = 1:swim_plot_stride:numel(x_swim_plot);
    green_swim = [0.47 0.76 0.36];
    pink_swim = [0.83 0.30 0.62];

    plot(x_swim_plot(swim_plot_idx), swim_base + swim_ch1(swim_plot_idx), ...
        'Color', green_swim, 'LineWidth', 0.85, ...
        'DisplayName', 'right swim');
    plot(x_swim_plot(swim_plot_idx), swim_base + swim_ch2(swim_plot_idx), ...
        'Color', pink_swim, 'LineWidth', 0.85, ...
        'DisplayName', 'left swim');
end

if ~isempty(stim_events)
    for event_id = 1:size(stim_events,1)
        onset_s = stim_events(event_id,1);
        offset_s = stim_events(event_id,2);
        dir_id = stim_events(event_id,3);
        mid_s = (onset_s + offset_s) / 2;
        theta = -pi/2 - (dir_id - 1) * pi / 6;

        line_x = [mid_s, mid_s + cos(theta) * direction_line_dx];
        line_y = [direction_line_y, direction_line_y + sin(theta) * direction_line_dy];
        plot(line_x, line_y, 'Color', trace_colors(dir_id,:), 'LineWidth', 2.3);
    end
end

ylim([direction_line_y - 0.14, 0.25]);
xlim(x_window);

scale_time_s = 10;
scale_dff = 0.10;
scale_x = x_window(2) - 50;
scale_y = -0.05;
scale_font_size = 25;
scale_text_gap = 5;
plot([scale_x - scale_time_s, scale_x], [scale_y, scale_y], 'k', 'LineWidth', 1.8);
plot([scale_x, scale_x], [scale_y, scale_y + scale_dff], 'k', 'LineWidth', 1.8);
text(scale_x + scale_text_gap, scale_y + scale_dff * 2, '10 %', ...
    'Color', 'k', 'FontSize', scale_font_size, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
text(scale_x + scale_text_gap, scale_y + scale_dff * 0.50, 'dff', ...
    'Color', 'k', 'FontSize', scale_font_size, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'Interpreter', 'tex');
text(scale_x - scale_time_s / 2, scale_y - 0.040, '10 s', ...
    'Color', 'k', 'FontSize', scale_font_size, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

set(gca, 'Visible', 'off', 'XTick', [], 'YTick', [], 'Box', 'off');
xlabel('');
ylabel('');

%%
route = strcat('D:\WSJ\Mulab\Paper_inbox\astroglia_direction\colorful_calcium_traces\',data_name,'_fig5.pdf');
route = Name_File_with_Suffix_astrodir2026(route);
exportgraphics(gcf, route, 'ContentType', 'vector', 'Resolution', 300);

%% plot direction-specific stimulus kernels
kernel_fig = figure('Color', 'w', 'Position', [100 100 980 420], 'Visible', 'off');
kernel_ax = axes(kernel_fig);
hold(kernel_ax, 'on');
set(kernel_ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
    'Box', 'off', 'TickDir', 'out');

kernel_height = 0.10;
frame_samples = round(x_scale_frames * SamplingFreq);

for m = 1:num_directions
    Mode1_ = [];
    Mode1_(2,:) = stages(m).onset;
    Mode1_(3,:) = stages(m).offset;
    Mode1_(1,:) = 1:size(Mode1_,2);

    Mode1_square_kernel = make_full_square_kernel_astrodir2026(Mode1_, Exp);
    kernel_trace = nan(size(x_scale_frames));
    valid_frame_mask = frame_samples >= 1 & frame_samples <= numel(Mode1_square_kernel);
    kernel_trace(valid_frame_mask) = Mode1_square_kernel(frame_samples(valid_frame_mask));
    kernel_trace = (kernel_trace - 1) / 0.05 * kernel_height;
    kernel_trace = kernel_trace - m * trace_spacing;

    plot(kernel_ax, x_scale_frames, kernel_trace, ...
        'Color', trace_colors(m,:), ...
        'LineWidth', 1.6);
end

xlim(kernel_ax, x_window);
ylim(kernel_ax, [-trace_spacing * (num_directions + 0.8), 0.05]);
set(kernel_ax, 'Visible', 'off', 'XTick', [], 'YTick', [], 'Box', 'off');
xlabel(kernel_ax, '');
ylabel(kernel_ax, '');

kernel_route = strcat('D:\WSJ\Mulab\Paper_inbox\astroglia_direction\colorful_calcium_traces\', ...
    data_name, '_kernels.pdf');
kernel_route = Name_File_with_Suffix_astrodir2026(kernel_route);
exportgraphics(kernel_fig, kernel_route, 'ContentType', 'vector', 'Resolution', 300);
close(kernel_fig);

function y = normalize_swim_trace_01(x, baseline_prctile, scale_prctile)
x = double(x(:));
finite_mask = isfinite(x);
if ~any(finite_mask)
    y = zeros(size(x));
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
