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

required_motion_files = {'ave.tif', 'motion_param.mat'};
missing_motion_files = required_motion_files(~cellfun(@(fn) exist(fullfile(data_dir, fn), 'file') == 2, required_motion_files));
if use_nb_data_release && ~isempty(missing_motion_files)
    warning('Release imaging directory is missing required motion file(s): %s. Falling back to original motion directory: %s', ...
        strjoin(missing_motion_files, ', '), legacy_data_dir);
    data_dir = legacy_data_dir;
end

plane_id = 12;
xy_pixel_size_um = 0.406;
scale_bar_um = 50;
scale_bar_px = scale_bar_um / xy_pixel_size_um;

patch_box_px = 20;
arrow_display_scale = 10;
show_scale_text = false;

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

ave_path = fullfile(data_dir, 'ave.tif');
motion_param_path = fullfile(data_dir, 'motion_param.mat');

load(motion_param_path, 'motion_param');

img = imread(ave_path, plane_id);
[img_h, img_w] = size(img);

mp = motion_param(plane_id);
inds = double(mp.indslist(:));
tilt_med = double(mp.tilt_med);

[patch_y, patch_x] = ind2sub([img_h, img_w], inds);

z_motion_rounded = round(tilt_med(:, 3));
z_values = [-2 -1 0 1 2];
z_colors = [
    0 0 1;
    0 1 1;
    0 1 0;
    1 1 0;
    1 0 0
];

z_pdf = fullfile(out_dir, sprintf('plane%02d_z_motion_patches.pdf', plane_id));
z_pdf = Name_File_with_Suffix_astrodir2026(z_pdf);
xy_pdf = fullfile(out_dir, sprintf('plane%02d_xy_motion_arrows.pdf', plane_id));
xy_pdf = Name_File_with_Suffix_astrodir2026(xy_pdf);

fig_z = new_full_image_figure(img_w, img_h);
ax_z = axes(fig_z, 'Position', [0 0 1 1]);
show_background(ax_z, img);
hold(ax_z, 'on');

half_box = patch_box_px / 2;
for i = 1:numel(z_values)
    idx = find(z_motion_rounded == z_values(i));
    for j = 1:numel(idx)
        k = idx(j);
        rectangle(ax_z, ...
            'Position', [patch_x(k) - half_box, patch_y(k) - half_box, patch_box_px, patch_box_px], ...
            'FaceColor', z_colors(i, :), ...
            'EdgeColor', [0 0 0], ...
            'LineWidth', 0.25);
    end
end

draw_scale_bar(ax_z, img_w, img_h, scale_bar_px, scale_bar_um, show_scale_text);
export_pdf(fig_z, z_pdf);
close(fig_z);

fig_xy = new_full_image_figure(img_w, img_h);
ax_xy = axes(fig_xy, 'Position', [0 0 1 1]);
show_background(ax_xy, img);
hold(ax_xy, 'on');

quiver(ax_xy, ...
    patch_x, patch_y, ...
    arrow_display_scale * tilt_med(:, 2), ...
    arrow_display_scale * tilt_med(:, 1), ...
    0, ...
    'Color', [1 0 0], ...
    'LineWidth', 1.0, ...
    'MaxHeadSize', 0.8);

draw_scale_bar(ax_xy, img_w, img_h, scale_bar_px, scale_bar_um, show_scale_text);
export_pdf(fig_xy, xy_pdf);
close(fig_xy);

fprintf('Saved z-motion PDF:  %s\n', z_pdf);
fprintf('Saved xy-motion PDF: %s\n', xy_pdf);
fprintf('Scale bar: %.1f um = %.2f pixels at %.3f um/pixel\n', ...
    scale_bar_um, scale_bar_px, xy_pixel_size_um);

function fig = new_full_image_figure(img_w, img_h)
fig = figure('Color', 'k', ...
    'Units', 'pixels', ...
    'Position', [100 100 img_w img_h], ...
    'InvertHardcopy', 'off', ...
    'Visible', 'off');
end

function show_background(ax, img)
clim = robust_clim(img, [0.5 100]);
imagesc(ax, img, clim);
colormap(ax, gray(256));
axis(ax, 'image');
axis(ax, 'off');
set(ax, 'YDir', 'reverse');
end

function draw_scale_bar(ax, img_w, img_h, scale_bar_px, scale_bar_um, show_scale_text)
margin_x = 70;
margin_y = 70;
bar_line_width = 7;

x2 = img_w - margin_x;
x1 = x2 - scale_bar_px;
y = margin_y;
y = img_h - margin_y;

plot(ax, [x1 x2], [y y], 'w-', 'LineWidth', bar_line_width, 'Clipping', 'off');

if show_scale_text
    text(ax, (x1 + x2) / 2, y - 28, sprintf('%g um', scale_bar_um), ...
        'Color', 'w', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 18, ...
        'FontWeight', 'bold');
end
end

function clim = robust_clim(img, pct)
x = double(img(:));
x = x(isfinite(x));
x = sort(x);
if isempty(x)
    clim = [0 1];
    return;
end

n = numel(x);
lo_idx = max(1, min(n, round(pct(1) / 100 * n)));
hi_idx = max(1, min(n, round(pct(2) / 100 * n)));
clim = [x(lo_idx), x(hi_idx)];

if clim(1) >= clim(2)
    clim = [min(x), max(x)];
end
if clim(1) >= clim(2)
    clim = [0 1];
end
end

function export_pdf(fig, pdf_path)
try
    exportgraphics(fig, pdf_path, ...
        'ContentType', 'vector', ...
        'BackgroundColor', 'black');
catch
    set(fig, 'PaperPositionMode', 'auto');
    print(fig, pdf_path, '-dpdf', '-painters');
end
end
