clear; close all; clc;

% Plot plane-12 Z- and XY-motion timecourses from the saved move_tcourse
% intermediate.
% The x-axis is converted from frame index to seconds using frames.mat.

data_dir = 'D:\WSJ\Mulab\Paper_inbox\astroglia_direction\motion_check_figures';
move_tcourse_path = fullfile(data_dir, 'move_tcourse_for_plot.mat');

use_nb_data_release = true;

if use_nb_data_release
    release_dir = 'F:\Data\Lightsheet\Astrocytes_direction\NB_data_release';
    data_name = 'fish_2'; % original data_name: 20240120_3_1
    frames_path = fullfile(release_dir, data_name, 'swimming', [data_name, '_processed'], [data_name, '_frames.mat']);
else
    frames_path = fullfile(data_dir, '20240120_3_1_frames.mat');
end

if 0
    use_sem = true;
    file_suffix = '_sem';
else
    use_sem = false;
    file_suffix = '_std';
end

z_pdf = fullfile(data_dir, ['motion_tcourse_z_shaded_seconds', file_suffix, '.pdf']);
xy_pdf = fullfile(data_dir, ['motion_tcourse_xy_shaded_seconds', file_suffix, '.pdf']);
z_pdf = Name_File_with_Suffix(z_pdf);
xy_pdf = Name_File_with_Suffix(xy_pdf);

plane_id = 12;
line_color = hex2rgb('#d5b34f');
shade_alpha = 0.28;
y_limits = [-10 10];

S = load(move_tcourse_path, 'move_tcourse', 'zpixeldist', 'xypixeldist');
F = load(frames_path, 'frames');

frame_interval_sec = (F.frames(2, 2) - F.frames(1, 2)) / 6000;

tcourse = double(S.move_tcourse(plane_id).tcourse(:));
time_sec = (tcourse - tcourse(1)) * frame_interval_sec;

z_mean_um = double(S.move_tcourse(plane_id).rs_ave_z(:)) * S.zpixeldist;
z_sd_um = double(S.move_tcourse(plane_id).rs_std_z(:)) * S.zpixeldist;

xy_mean_um = double(S.move_tcourse(plane_id).rs_ave_xy(:)) * S.xypixeldist;
xy_sd_um = double(S.move_tcourse(plane_id).rs_std_xy(:)) * S.xypixeldist;

n_patches = double(S.move_tcourse(plane_id).n_patches);
if use_sem
    z_err_um = z_sd_um / sqrt(n_patches);
    xy_err_um = xy_sd_um / sqrt(n_patches);
    error_label = 'SEM';
else
    z_err_um = z_sd_um;
    xy_err_um = xy_sd_um;
    error_label = 'SD';
end

plot_shaded_timecourse(time_sec, z_mean_um, z_err_um, ...
    'Z displacement (\mum)', y_limits, line_color, shade_alpha, z_pdf);

plot_shaded_timecourse(time_sec, xy_mean_um, xy_err_um, ...
    'XY displacement (\mum)', y_limits, line_color, shade_alpha, xy_pdf);

fprintf('Saved plane-%02d Z shaded timecourse PDF: %s\n', plane_id, z_pdf);
fprintf('Saved plane-%02d XY shaded timecourse PDF: %s\n', plane_id, xy_pdf);
fprintf('Shaded error: %s, n_patches=%d\n', error_label, n_patches);
fprintf('Frame interval: %.6f s\n', frame_interval_sec);


function rgb = hex2rgb(hex)
hex = char(hex);
if startsWith(hex, '#')
    hex = hex(2:end);
end
rgb = sscanf(hex, '%2x%2x%2x', [1 3]) / 255;
end


function plot_shaded_timecourse(time_sec, mean_um, err_um, y_label, y_limits, ...
    line_color, shade_alpha, out_pdf)
upper = mean_um + err_um;
lower = mean_um - err_um;

fig = figure('Color', 'w', ...
    'Units', 'pixels', ...
    'Position', [250 250 620 460], ...
    'Visible', 'off');
ax = axes(fig, 'Position', [0.16 0.16 0.78 0.76]);
hold(ax, 'on');

fill(ax, ...
    [time_sec; flipud(time_sec)], ...
    [upper; flipud(lower)], ...
    line_color, ...
    'FaceAlpha', shade_alpha, ...
    'EdgeColor', 'none');

plot(ax, time_sec, mean_um, ...
    'Color', line_color, ...
    'LineWidth', 2.0);

xlim(ax, [0 max(time_sec)]);
ylim(ax, y_limits);
% ylim(ax, [-abs(min(upper)) * 4, max(upper) * 4]);
xlabel(ax, 'Time (s)');
ylabel(ax, y_label);

style_axes(ax);

exportgraphics(fig, out_pdf, ...
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

grid(ax, 'off');
ax.XAxis.FontSize = 11;
ax.YAxis.FontSize = 11;
ax.XLabel.FontSize = 13;
ax.YLabel.FontSize = 13;
end
