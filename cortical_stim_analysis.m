clear all
close all
clc

%% Load Data
folder_name = string('G:\.shortcut-targets-by-id\1DLmzFkqnKdtnmtAAuam9z99HmOFaf-OA\scs_mapping\scgs004\scgs004_2022-01-10_emg-implant\ephys\scgs004_2021-01-10_05_cxs_mur')
split_name = split(folder_name,'\')
size_name = numel(split_name);
new_name = '';
for i=1:numel(split_name);
    if i==size_name;
        break
    end
    if i==1;
        new_name = strcat (new_name, split_name(i))
    else
        new_name = strcat (new_name,'\', split_name(i))
    end
end

cd (folder_name)
addpath(genpath(' C:\Users\asa2248\Desktop\mapping\TDTMatlabSDK (1)\TDTSDK\TDTbin2mat')); %Add with subfolders
addpath(genpath(new_name))
data = TDTbin2mat(char(folder_name))

emg_data = data.streams.MUw4.data; 
biceps = double(emg_data(1,:));
ecr = double(emg_data(2,:));
fcr = double(emg_data(3,:));
biceps_right = double(emg_data(8,:));

%% Filter data
fl = 130;
fu = 1000;
fs = data.streams.MUw4.fs;

[b,a] = butter(2, [(2*fl)/fs (2*fu)/fs],'bandpass');
final_biceps = filtfilt(b,a,biceps);
final_ecr = filtfilt(b,a,ecr);
final_fcr = filtfilt(b,a,fcr);
final_biceps_right = filtfilt(b,a,biceps_right);


%% extract the stimulation/response window
% cd 'G:\.shortcut-targets-by-id\1DLmzFkqnKdtnmtAAuam9z99HmOFaf-OA\scs_mapping\scgs002\scgs002_2021-11-30_scs-implant\ephys\scgs002_2021-11-30_01_scsmap'
stim_time = data.scalars.eS1p.ts;
direct = dir('*.mat');
stim_table = load (direct.name);
stim_num = size(stim_table.T);
num_of_stim = stim_num(1); % Opens table to check the total number of stim

stim_intensities = unique(stim_table.T.pulse_amplitude);
count_int = histc(stim_table.T.pulse_amplitude,stim_intensities);

snap_biceps = {};
snap_ecr = {};
snap_fcr = {};
snap_biceps_right = {};
window_len = 0.05;
for i=1:num_of_stim; 
    
    snap_biceps = [snap_biceps [final_biceps((stim_time(i)-window_len)*fs : (stim_time(i)+window_len)*fs)]];
    snap_ecr = [snap_ecr [final_ecr((stim_time(i)-window_len)*fs : (stim_time(i)+window_len)*fs)]];
%     snap_ecr{end+1} = [final_ecr((stim_time(i)-window_len)*fs : (stim_time(i)+window_len)*fs)];
    snap_fcr = [snap_fcr [final_fcr((stim_time(i)-window_len)*fs : (stim_time(i)+window_len)*fs)]];
    snap_biceps_right = [snap_biceps_right [final_biceps_right((stim_time(i)-window_len)*fs : (stim_time(i)+window_len)*fs)]];

end

%% plot time
tt = numel(snap_biceps{1})/fs;
time = -tt/2:1/fs:(tt/2-(1/fs));

%% average and plot
reshaped_snap_biceps = reshape(snap_biceps,[numel(snap_biceps)/count_int(1),count_int(1)]);
reshaped_snap_ecr = reshape(snap_ecr,[numel(snap_ecr)/count_int(1),count_int(1)]);
reshaped_snap_fcr = reshape(snap_fcr,[numel(snap_fcr)/count_int(1),count_int(1)]);
reshaped_snap_biceps_right = reshape(snap_biceps_right,[numel(snap_biceps_right)/count_int(1),count_int(1)]);

snap_biceps_ave = zeros(1,numel(reshaped_snap_biceps{1}));
snap_ecr_ave = zeros(1,numel(reshaped_snap_ecr{1}));
snap_fcr_ave = zeros(1,numel(reshaped_snap_fcr{1}));
snap_right_biceps_ave = zeros(1,numel(reshaped_snap_biceps_right{1}));

for i=1:length(reshaped_snap_biceps)
    for ii = 1:count_int(1);
        snap_biceps_ave = snap_biceps_ave+reshaped_snap_biceps{i,ii}; 
        snap_ecr_ave = snap_ecr_ave+reshaped_snap_ecr{i,ii}; 
        snap_fcr_ave = snap_fcr_ave+reshaped_snap_fcr{i,ii}; 
        snap_right_biceps_ave = snap_right_biceps_ave+reshaped_snap_biceps_right{i,ii}; 
    end
    snap_biceps_final{i}=snap_biceps_ave/count_int(1);
    snap_ecr_final{i}=snap_ecr_ave/count_int(1);
    snap_fcr_final{i}=snap_fcr_ave/count_int(1);
    snap_biceps_right_final{i}=snap_right_biceps_ave/count_int(1);
end

figure
subplot(1,4,1)
b=0;
b_inc = 0.0001;
for i=1:numel(snap_biceps_final)
    subplot(1,4,1);plot(time,snap_biceps_final{i}+b);hold on
    subplot(1,4,2);plot(time,snap_ecr_final{i}+b);hold on
    subplot(1,4,3);plot(time,snap_fcr_final{i}+b);hold on
    subplot(1,4,4);plot(time,snap_biceps_right_final{i}+b);hold on
    b=b+b_inc
end

upper_lim = 20e-4;
lower_lim = 2e-4;
subplot(1,4,1);box off;title('Biceps');ylabel('Amplitude (uA)');ylim([-lower_lim upper_lim ])
yticklabels(num2cell(stim_intensities))
yticks([0:b_inc:(numel(stim_intensities)-1)*b_inc])
subplot(1,4,2);box off;title('ECR');yticklabels([]);ylim([-lower_lim upper_lim ])
subplot(1,4,3);box off;title('FCR');yticklabels([]);ylim([-lower_lim upper_lim ])
subplot(1,4,4);box off;title('Biceps right');yticklabels([]);ylim([-lower_lim upper_lim ])



