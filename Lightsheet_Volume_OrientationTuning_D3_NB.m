% pipeline for analyzing the astrocytic response to directional visual-motion stimulus dataset

%% import ephys data
disp('Clearing workspace.');
clear all
close all

tic
SamplingFreq = 6000;
data_name = '20240121_2_1_6d_gfap_g6f_D3_orientation';
date = string(data_name(:,1:8));
data_id = string(data_name(:,1:12));

data_root = 'F:\Data\Lightsheet\Astrocytes_direction\';
VR_dir = convertStringsToChars(strcat(data_root,date,'\VR\'));
Image_dir = convertStringsToChars(strcat(data_root,date,'\Volume\'));
data = read_data_WSJ_D2(strcat(VR_dir,data_name,'.10chFlt'));


x_scale = (1:length(data.ch2))/6000;
filtdata1 = filter_data(data.ch1);
filtdata2 = filter_data(data.ch2);

figure;hold on;plot(x_scale(1:50:end),filtdata1(1:50:end),'r','DisplayName','Channel 1')
plot(x_scale(1:50:end),-filtdata2(1:50:end),'b','DisplayName','Channel 2')
legend
ylim([-0.005 0.005]);

% choose which channel to filter
[~,y,button] = ginput(1);
if y > 0 % Double-press the space bar to keep the current trial; double-click the left mouse button to discard it.
    filtdata = filtdata1;
    disp("Using channel 1 ephys data.");
elseif y < 0
    filtdata = filtdata2;
    disp("Using channel 2 ephys data.");
end

raw_data = data;
raw_filtdata = filtdata;
figure;hold on;plot(x_scale(1:50:end),filtdata(1:50:end));
plot(x_scale(1:200:end),data.stimGain(1:200:end)/20,'color','#DF3A3E','linewidth',1.5);
toc

%% plot results
close all
figure;hold on
plot(x_scale(1:50:end),filtdata(1:50:end)*200);
% h1 = plot(x_scale(1:200:end),data.Trial_Mode(1:200:end)); %/100,'r','linewidth',1.5);
% plot(x_scale(1:100:end),data.fish_vel(1:100:end)*150);
plot(x_scale(1:100:end), data.orient(1:100:end)/100);
% plot(x_scale(1:200:end),data.stimGain(1:200:end)*10);
% plot(x_scale(1:200:end),data.Trial_Stage(1:200:end)/100);
% set(gca,'fontsize',50,'fontname','Arial','FontWeight','bold','Color','none');
% set(gca,'XColor','w');
% set(gca,'YColor','w');
% set(gcf,'color','none');
ylim([0 0.5]);

%% load cell response data
tic
close all
disp("Loading imaging data...")
image_dir = convertStringsToChars(strcat(Image_dir,string(data_id),'_registered'));
cell_resp_dim_processed = [];
load(fullfile(image_dir,'\cell_resp_dim_processed.mat'));
cell_resp_processed = [];
cell_resp_processed = read_LSstack_fast_float(fullfile(image_dir,'\cell_resp_processed.stackf'),cell_resp_dim);
disp('cell_resp_processed loaded');
cell_resp_raw = [];
cell_resp_raw = read_LSstack_fast_float(fullfile(image_dir,'\cell_resp.stackf'),cell_resp_dim);
disp('cell_resp_raw loaded');
cell_resp_lowcut = [];
cell_resp_lowcut = read_LSstack_fast_float(fullfile(image_dir,'\cell_resp_lowcut.stackf'),cell_resp_dim);
disp('cell_resp_lowcut loaded');
cell_info = [];
if exist(fullfile(image_dir,'\cell_loc.csv'))
    load(fullfile(image_dir,'\cell_info_processed.mat'));
    disp('Registered.');
    cells = []; 
    for k = 1: cell_resp_dim(1)
        cells(k,1) = cell_info(k).center_aligned(2); % If this errors, cell2region.m may not have been run.
        cells(k,2) = cell_info(k).center_aligned(1);
        cells(k,3) = cell_info(k).center_aligned(3);
    end
else
    load(fullfile(image_dir,'\cell_info.mat'));
    disp('Not registered.'); 
end

% load planes and exposure time
xml_path = strcat(image_dir,'\ch0_cam1.xml');
if contains(xml_path, 'cam1')
    cam_id = 2;
elseif contains(xml_path, 'cam0')
    cam_id = 1;
end
camera_settings = Extract_Parameters_xml(xml_path,'exposure_ms','dimensions');
exposure_time = camera_settings{1};
planes = cell2mat(camera_settings{1,2});
planes = planes(end);
m = 1; % First frame to read
n = cell_resp_dim(2);
% n = 6600;

toc

%% extract frames between imaging
tic
frames = []; % Store C# sample indices for each frame: column 1 is frame ID, column 2 is the C# sample index.

% for volumetric
if contains(xml_path, 'cam0') % jrgeco
    frames = find(data.ch4(1:end-1) < 3.6 & data.ch4(2:end) > 3.6); % Find the ID of each stack.
elseif contains(xml_path, 'cam1') % gcamp
    frames = find(data.ch3(1:end-1) < 3.6 & data.ch3(2:end) > 3.6); % Find the ID of each stack.
end

imaging_freq_by_stack = 1/((frames(2) - frames(1))/6000); % Estimate stack sampling frequency for rough plotting only; do not use it to identify stage onset frames.

% frames = frames(1:1599); % 20230325_1_1
frame_start = frames(1);
frame_end = frames(length(frames)); % frames is still a one-dimensional vector here.

if exist('m','var') && ~exist('n','var')
    frame_start = frames(m);
    frames(1:(m-1)) = [];
elseif exist('n','var') && ~exist('m','var')
    frame_end = frames(n);
    frames((n+1):length(frames),:) = [];
elseif exist('n','var') && exist('m','var')
    frame_start = frames(m);
    frame_end = frames(n);
    frames(1:(m-1)) = 0;
    frames((n+1):length(frames),:) = 0;
    frames(all(~frames,2),:) = []; 
end

frames(:,2) = frames(:,1); % C# sample index
frames(:,1) = 1:size(frames,1); % frame id

if size(frames,1) == size(cell_resp_raw,2) + 1
    frames(size(frames,1),:) = [];
    frame_end = frames(size(frames,1),2);
    disp("Frames aligned.")
end

frames(:,2) = frames(:,2) - frame_start + 1; 
x_scale_frames = [];
x_scale_frames = (frames(:,2)/6000).';

filtdata = [];
filtdata = raw_filtdata(frame_start:frame_end);
x_scale = [];
x_scale = (1:length(filtdata))/6000;

% align each fields of this structure based on frame start
fields = [];
fields = fieldnames(data);
val = [];
for i = 1:numel(fields)
   val = data.(fields{i});
   data.(fields{i}) = val(frame_start:frame_end);
end

StimGain_frames = [];
StimGain_frames = data.stimGain(frames(:,2)).'; % get the stim gain in each frame
x_scale_stimgain = [];
x_scale_stimgain = (frames(:,2)/6000).';

toc

%% extract online swims and offline swims
tic
Bouts_offline = {};
Bouts_offline = extract_bouts(filtdata);

% distinguish OL swim and CL swim in each trial
% by detecting whether velcosity breaks zero line within each bout
Bout_Types = [];
Bout_Types = find_swim_type(Bouts_offline,data.fish_vel,min(Bouts_offline.bout_end - Bouts_offline.bout_start));

% assign bout type to each bout
Bouts_offline.bout_type = zeros(1,length(Bouts_offline.bout_start));
Bouts_offline.bout_type(Bout_Types.CL_swim) = 1;

% check bout types
figure;hold on
xlim([100 500]);
ylim([-0.001 0.01])
p1 = plot(x_scale(1:100:end),filtdata(1:100:end));
p2 = plot((Bouts_offline.bout_start(Bout_Types.CL_swim))/6000,0.001,'r*','DisplayName','CL swim');
p3 = plot((Bouts_offline.bout_start(Bout_Types.OL_swim))/6000,0.00001,'b*','DisplayName','OL swim');
p4 = plot(x_scale(1:100:end),data.fish_vel(1:100:end),'color','#F6821F');
toc

