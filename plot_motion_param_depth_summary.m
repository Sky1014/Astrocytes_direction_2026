script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'astro_functions'));

clear; close all; clc;

use_nb_data_release = true;
legacy_data_dir = 'F:\Data\Lightsheet\20240120\Volume\20240120_3_1_registered';

if use_nb_data_release
    release_dir = 'F:\Data\Lightsheet\Astrocytes_direction\NB_data_release';
    data_name = 'fish_2';
    data_dir = fullfile(release_dir, data_name, 'imaging', [data_name, '_registered']);
else
    data_dir = legacy_data_dir;
end
out_dir = 'D:\WSJ\Mulab\Paper_inbox\astroglia_direction\motion_check_figures';

required_motion_files = {'motion_param.mat'};
missing_motion_files = required_motion_files(~cellfun(@(fn) exist(fullfile(data_dir, fn), 'file') == 2, required_motion_files));
if use_nb_data_release && ~isempty(missing_motion_files)
    warning('Release imaging directory is missing required motion file(s): %s. Falling back to original motion directory: %s', ...
        strjoin(missing_motion_files, ', '), legacy_data_dir);
    data_dir = legacy_data_dir;
end

motion_param_path = fullfile(data_dir, 'motion_param.mat');
xy_pdf = fullfile(out_dir, 'motion_param_xy_depth_summary.pdf');
z_pdf = fullfile(out_dir, 'motion_param_z_depth_summary.pdf');
xy_pdf = Name_File_with_Suffix_astrodir2026(xy_pdf);
z_pdf = Name_File_with_Suffix_astrodir2026(z_pdf);

depth_step_um = 8;
xy_pixel_size_um = 0.406;
z_step_um = 5;

blue = [0.35 0.68 0.88];
blue_edge = [0.10 0.44 0.67];

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

load(motion_param_path, 'motion_param');

n_planes = numel(motion_param);
depth_um = (0:n_planes-1) * depth_step_um;

xy_mean = nan(1, n_planes);
xy_sem = nan(1, n_planes);
z_mean = nan(1, n_planes);
z_sem = nan(1, n_planes);

for zz = 1:n_planes
    tilt_med = double(motion_param(zz).tilt_med);

    xy_values = sqrt(tilt_med(:, 1).^2 + tilt_med(:, 2).^2) * xy_pixel_size_um;
    z_values = abs(tilt_med(:, 3)) * z_step_um;

    xy_mean(zz) = local_mean_omitnan(xy_values);
    xy_sem(zz) = local_sem_omitnan(xy_values);
    z_mean(zz) = local_mean_omitnan(z_values);
    z_sem(zz) = local_sem_omitnan(z_values);
end

plot_depth_summary(depth_um, xy_mean, xy_sem, ...
    'Horizontal displacement (\mum)', 'XY motion summary', ...
    xy_pdf, blue, blue_edge, [-1 5]);

plot_depth_summary(depth_um, z_mean, z_sem, ...
    'Axial displacement (\mum)', 'Z motion summary', ...
    z_pdf, blue, blue_edge, [-1 10]);

fprintf('Saved XY summary PDF: %s\n', xy_pdf);
fprintf('Saved Z summary PDF:  %s\n', z_pdf);

function plot_depth_summary(depth_um, mean_values, sem_values, ...
    y_label, fig_title, pdf_path, blue, blue_edge, y_limits)

fig = figure('Color', 'w', ...
    'Units', 'pixels', ...
    'Position', [250 250 520 620], ...
    'Visible', 'off');
ax = axes(fig, 'Position', [0.16 0.14 0.78 0.80]);
hold(ax, 'on');

bar_relative_width = 0.65;
bar(ax, depth_um, mean_values, ...
    bar_relative_width, ...
    'FaceColor', blue, ...
    'EdgeColor', blue_edge, ...
    'LineWidth', 1.15, ...
    'FaceAlpha', 0.70);

errorbar(ax, depth_um, mean_values, sem_values, ...
    'LineStyle', 'none', ...
    'Color', blue_edge, ...
    'LineWidth', 1.6, ...
    'CapSize', 8);

x_axis_max = max(depth_um) + depth_step_from_values(depth_um);
xlim(ax, [-2, x_axis_max]);
ylim(ax, y_limits);

xticks(ax, 0:16:x_axis_max);
xlabel(ax, 'Imaging depth (\mum)');
ylabel(ax, y_label);

style_axes(ax);

exportgraphics(fig, pdf_path, ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');
close(fig);
end

function style_axes(ax)
set(ax, ...
    'Box', 'on', ...
    'TickDir', 'out', ...
    'TickLength', [0.012 0.012], ...
    'LineWidth', 0.75, ...
    'FontName', 'Arial', ...
    'FontSize', 12, ...
    'XColor', [0.25 0.25 0.25], ...
    'YColor', [0.25 0.25 0.25]);

ax.XAxis.FontSize = 11;
ax.YAxis.FontSize = 11;
ax.XLabel.FontSize = 13;
ax.YLabel.FontSize = 13;
ax.Title.FontSize = 12;
end

function y = local_mean_omitnan(x)
x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    y = mean(x);
end
end

function y = local_sem_omitnan(x)
x = x(isfinite(x));
n = numel(x);
if n <= 1
    y = NaN;
else
    y = std(x, 0) / sqrt(n);
end
end

function step = depth_step_from_values(depth_um)
if numel(depth_um) < 2
    step = 8;
else
    step = median(diff(depth_um));
end
end
