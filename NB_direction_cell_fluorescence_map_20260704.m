script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'astro_functions'));

disp('Clearing workspace.');
clear;
close all;

%% Paths and user options

preference_method = 'regressor';

num_directions = 12;
rho_standard = 0.1;
min_peak_dff = 0.02;
marker_size = 8;
marker_alpha = 0.88;
background_brightness = 0.55;
show_direction_key = false;
show_title = false;
save_figure = false;

export_route = 'D:\WSJ\Mulab\Paper_inbox\astroglia_direction\NB_direction_cell_map.pdf';
export_route = Name_File_with_Suffix_astrodir2026(export_route);

%% Import data, following D3_NB_Figure_withSwim.m
tic;
use_nb_data_release = true;

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

    VR_dir = convertStringsToChars(strcat(raw_dir, date, '\VR\'));
    Image_dir = convertStringsToChars(strcat(raw_dir, date, '\Volume\'));
    folder = strcat(VR_dir, data_id, '_processed');
    image_dir = convertStringsToChars(strcat(Image_dir, string(data_id), '_registered'));
end

if use_nb_data_release
    load(strcat(folder, '\', data_id, '_right_swim.mat'), 'right_swim');
    load(strcat(folder, '\', data_id, '_left_swim.mat'), 'left_swim');
    filtdata1 = right_swim;
    filtdata2 = left_swim;
else
    load(strcat(folder, '\', data_id, '_ch1.mat'));
    load(strcat(folder, '\', data_id, '_ch2.mat'));
    filtdata1 = filter_data_astrodir2026(ch1);
    filtdata2 = filter_data_astrodir2026(ch2);
end

load(strcat(folder, '\', data_id, '_orient.mat'));
load(strcat(folder, '\', data_id, '_frames.mat'));
load(strcat(folder, '\', data_id, '_trials.mat'));
load(strcat(folder, '\', data_id, '_stages.mat'));
load(strcat(folder, '\', data_id, '_stages_framed.mat'));

disp('Loading imaging data...');
cell_resp_dim_processed = [];
load(fullfile(image_dir, 'cell_resp_dim_processed.mat'));
cell_resp_processed = read_LSstack_fast_float_astrodir2026(fullfile(image_dir, 'cell_resp_processed.stackf'), cell_resp_dim);
disp('cell_resp_processed loaded.');

cell_info = [];
registered_cells = exist(fullfile(image_dir, 'cell_loc.csv'), 'file') == 2;

if registered_cells
    load(fullfile(image_dir, 'cell_info_processed.mat'));
    disp('Registered cell locations loaded.');
    cells = zeros(cell_resp_dim(1), 3);
    for k = 1:cell_resp_dim(1)
        cells(k,1) = cell_info(k).center_aligned(2);
        cells(k,2) = cell_info(k).center_aligned(1);
        cells(k,3) = cell_info(k).center_aligned(3);
    end
else
    load(fullfile(image_dir, 'cell_info.mat'));
    disp('Unregistered cell locations loaded.');
    cells = zeros(cell_resp_dim(1), 3);
    for k = 1:cell_resp_dim(1)
        cells(k,1) = cell_info(k).center(2);
        cells(k,2) = cell_info(k).center(1);
        cells(k,3) = cell_info(k).slice;
    end
end

%% Camera and background brain
cam_id = 2;

if registered_cells
    template_dir = 'D:\WSJ\Code\git\Multineuromodulatory-integration\13_registration\gfap\';
    template_path = [template_dir, '\Temp_gfapChR2ECFP_8bit.nrrd'];
    brain_bg = double(nrrdread_astrodir2026(template_path));
else
    ave_path = [image_dir, '\ave.tif'];
    brain_bg = double(ReadTiff_astrodir2026(ave_path));
end

%% Direction colors, following D3_NB_Figure_withSwim.m
tmp_fig = figure('Visible', 'off');
colormap_glia = colormap(slanCM_astrodir2026('gist_rainbow'));
close(tmp_fig);

row_ids = round(linspace(1, size(colormap_glia,1), num_directions));
row_ids(row_ids < 1) = 1;
row_ids(row_ids > size(colormap_glia,1)) = size(colormap_glia,1);
direction_colors = colormap_glia(row_ids, :);

direction_hsv = rgb2hsv(direction_colors);
direction_hsv(:,2) = direction_hsv(:,2) * 0.60;
direction_hsv(:,3) = direction_hsv(:,3) * 0.82;
direction_colors = hsv2rgb(direction_hsv);

%% Estimate preferred direction for each cell
valid_cells = all(isfinite(cell_resp_processed), 2);
valid_cell_ids = find(valid_cells);

switch lower(preference_method)
    case 'regressor'
        [preferred_dir, preference_score, direction_score] = local_regressor_preference( ...
            cell_resp_processed, valid_cell_ids, frames, stages, num_directions, rho_standard);
        selected_cells = find(isfinite(preferred_dir) & preference_score >= rho_standard);

    case 'trial_peak'
        [preferred_dir, preference_score, direction_score] = local_trial_peak_preference( ...
            cell_resp_processed, valid_cell_ids, stages_framed, num_directions, min_peak_dff);
        selected_cells = find(isfinite(preferred_dir) & preference_score >= min_peak_dff);

    otherwise
        error('Unknown preference_method: %s', preference_method);
end

fprintf('Selected %d / %d cells using %s preference.\n', ...
    numel(selected_cells), size(cell_resp_processed,1), preference_method);

NB_direction_map = struct();
NB_direction_map.data_name = data_name;
NB_direction_map.image_dir = image_dir;
NB_direction_map.preference_method = preference_method;
NB_direction_map.rho_standard = rho_standard;
NB_direction_map.min_peak_dff = min_peak_dff;
NB_direction_map.preferred_dir = preferred_dir;
NB_direction_map.preference_score = preference_score;
NB_direction_map.direction_score = direction_score;
NB_direction_map.selected_cells = selected_cells;
NB_direction_map.direction_colors = direction_colors;

%% Plot top-view fluorescence-style map
[bg_projection, plot_cells] = local_prepare_top_view(brain_bg, cells, cam_id);
bg_projection = local_contrast_scale(bg_projection, 1, 99.8) * background_brightness;

plot_cells = plot_cells(selected_cells, :);
plot_dir = preferred_dir(selected_cells);
plot_score = preference_score(selected_cells);

[plot_score, order] = sort(plot_score, 'ascend');
plot_cells = plot_cells(order, :);
plot_dir = plot_dir(order);
selected_cells = selected_cells(order);

cell_rgb = direction_colors(plot_dir, :);

fig = figure('Color', 'k', 'Position', [100 80 560 980]);
ax = axes(fig);
imagesc(ax, bg_projection, [0 1]);
colormap(ax, gray);
axis(ax, 'image');
axis(ax, 'ij');
axis(ax, 'off');
hold(ax, 'on');

scatter_handle = scatter(ax, plot_cells(:,1), plot_cells(:,2), marker_size, cell_rgb, ...
    'filled', 'MarkerEdgeColor', 'none');

try
    scatter_handle.MarkerFaceAlpha = marker_alpha;
catch
    warning('Marker alpha is not supported in this MATLAB version.');
end

set(ax, 'Color', 'k');
set(fig, 'InvertHardCopy', 'off');

if show_title
    title(ax, sprintf('NB direction map, %s, n = %d', preference_method, numel(selected_cells)), ...
        'Color', 'w', 'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'normal');
end

if show_direction_key
    local_plot_direction_key(ax, bg_projection, direction_colors);
end

if save_figure
    exportgraphics(fig, export_route, 'ContentType', 'image', 'Resolution', 600);
    disp(['Figure saved to: ', export_route]);
end

toc;

%% Local functions
function [preferred_dir, preference_score, direction_score] = local_regressor_preference( ...
    cell_resp_processed, valid_cell_ids, frames, stages, num_directions, rho_standard)

    direction_score = nan(size(cell_resp_processed,1), num_directions);

    Exp = [];
    Exp(1,1) = 1;
    Exp(2,1) = stages(1).onset(1);
    last_stage_id = min(11, numel(stages));
    Exp(3,1) = stages(last_stage_id).offset(end);

    frame_mask = frames(:,2) >= Exp(2,1) & frames(:,2) <= Exp(3,1);
    frame_qualified = frames(frame_mask, 2);
    cell_resp_qualified = cell_resp_processed(valid_cell_ids, frame_mask);
    frame_scale_qualified = (frame_qualified - frame_qualified(1) + 1) / 6000;

    for dir_id = 1:num_directions
        Mode = [];
        Mode(2,:) = stages(dir_id).onset;
        Mode(3,:) = stages(dir_id).offset;
        Mode(1,:) = 1:size(Mode,2);

        Mode_square_kernel = make_full_square_kernel_astrodir2026(Mode, Exp);
        Mode_cells = make_customized_regressor_astrodir2026(Mode_square_kernel, frame_scale_qualified, ...
            frame_qualified, cell_resp_qualified, rho_standard);

        original_cell_ids = valid_cell_ids(Mode_cells.cell_index);
        direction_score(original_cell_ids, dir_id) = Mode_cells.rho(:);
    end

    score_for_max = direction_score;
    score_for_max(~isfinite(score_for_max)) = -Inf;
    [preference_score, preferred_dir] = max(score_for_max, [], 2);
    preferred_dir(~isfinite(preference_score) | preference_score == -Inf) = NaN;
    preference_score(preference_score == -Inf) = NaN;
end

function [preferred_dir, preference_score, direction_score] = local_trial_peak_preference( ...
    cell_resp_processed, valid_cell_ids, stages_framed, num_directions, min_peak_dff)

    direction_score = nan(size(cell_resp_processed,1), num_directions);

    for dir_id = 1:num_directions
        if dir_id > numel(stages_framed) || isempty(stages_framed(dir_id).framed_onset)
            continue;
        end

        trial_peak = nan(numel(valid_cell_ids), numel(stages_framed(dir_id).framed_onset));
        for trial_id = 1:numel(stages_framed(dir_id).framed_onset)
            onset_frame = max(1, stages_framed(dir_id).framed_onset(trial_id));
            offset_frame = min(size(cell_resp_processed,2), stages_framed(dir_id).framed_offset(trial_id));
            if offset_frame <= onset_frame
                continue;
            end
            trial_peak(:, trial_id) = max(cell_resp_processed(valid_cell_ids, onset_frame:offset_frame) - 1, [], 2);
        end

        direction_score(valid_cell_ids, dir_id) = local_nanmean(trial_peak, 2);
    end

    score_for_max = direction_score;
    score_for_max(~isfinite(score_for_max)) = -Inf;
    [preference_score, preferred_dir] = max(score_for_max, [], 2);
    preferred_dir(~isfinite(preference_score) | preference_score == -Inf) = NaN;
    preference_score(preference_score == -Inf) = NaN;
    preferred_dir(preference_score < min_peak_dff) = NaN;
end

function [bg_projection, plot_cells] = local_prepare_top_view(brain_bg, cells, cam_id)
    bg_projection = squeeze(mean(brain_bg, 3));
    plot_cells = cells;

    if cam_id == 2
        plot_cells(:,1) = size(bg_projection, 2) - plot_cells(:,1);
    elseif cam_id == 1
        bg_projection = fliplr(bg_projection);
    else
        error('Unsupported cam_id: %d', cam_id);
    end
end

function scaled_img = local_contrast_scale(img, low_percentile, high_percentile)
    img = double(img);
    vals = img(isfinite(img));
    vals = sort(vals(:));

    if isempty(vals)
        scaled_img = zeros(size(img));
        return;
    end

    low_idx = max(1, round(low_percentile / 100 * numel(vals)));
    high_idx = min(numel(vals), round(high_percentile / 100 * numel(vals)));
    low_val = vals(low_idx);
    high_val = vals(high_idx);

    if high_val <= low_val
        high_val = max(vals);
        low_val = min(vals);
    end

    if high_val <= low_val
        scaled_img = zeros(size(img));
    else
        scaled_img = (img - low_val) / (high_val - low_val);
        scaled_img(scaled_img < 0) = 0;
        scaled_img(scaled_img > 1) = 1;
    end
end

function y = local_nanmean(x, dim)
    finite_mask = isfinite(x);
    x(~finite_mask) = 0;
    n = sum(finite_mask, dim);
    y = sum(x, dim) ./ n;
    y(n == 0) = NaN;
end

function local_plot_direction_key(ax, bg_projection, direction_colors)
    x0 = size(bg_projection, 2) * 0.08;
    y0 = size(bg_projection, 1) * 0.08;
    radius = min(size(bg_projection)) * 0.035;
    gap = radius * 2.2;

    for dir_id = 1:size(direction_colors,1)
        x = x0 + mod(dir_id - 1, 3) * gap;
        y = y0 + floor((dir_id - 1) / 3) * gap;
        scatter(ax, x, y, 35, direction_colors(dir_id,:), 'filled', ...
            'MarkerEdgeColor', 'none');
        text(ax, x + radius * 0.65, y, num2str(dir_id), 'Color', 'w', ...
            'FontName', 'Arial', 'FontSize', 8, 'VerticalAlignment', 'middle');
    end
end
