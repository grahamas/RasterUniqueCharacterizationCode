clearvars
close all
clc

%load data

datafiles = {
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_A2_trial1.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_A3_trial2.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_A4_trial3.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_A5_trial4.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_A6_trial5.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_A7_trial6.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_A8_trial7.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_B2_trial8.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_B3_trial9.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_B4_trial10.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_B5_trial11.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_B6_trial12.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_B7_trial13.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_B8_trial14.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_C2_trial15.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_C3_trial16.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_C4_trial17.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_C5_trial18.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_C6_trial19.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_C7_trial20.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_C8_trial21.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_D2_trial22.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_D3_trial23.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_D4_trial24.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_D5_trial25.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_D6_trial26.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_D7_trial27.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_D8_trial28.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_E2_trial29.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_E3_trial30.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_E4_trial31.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_E5_trial32.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_E6_trial33.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_E7_trial34.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_E8_trial35.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_F2_trial36.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_F3_trial37.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_F4_trial38.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_F5_trial39.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_F6_trial40.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_F7_trial41.mat',
    'Rat Baseline Dataset/Rat_MEA2Baseline_DIV22_raster_3D_F8_trial42.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_A2_trial1.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_A3_trial2.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_A4_trial3.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_A5_trial4.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_A6_trial5.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_A7_trial6.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_A8_trial7.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_B2_trial8.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_B3_trial9.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_B4_trial10.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_B5_trial11.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_B6_trial12.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_B7_trial13.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_B8_trial14.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_C2_trial15.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_C3_trial16.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_C4_trial17.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_C5_trial18.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_C6_trial19.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_C7_trial20.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_C8_trial21.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_D2_trial22.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_D3_trial23.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_D4_trial24.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_D5_trial25.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_D6_trial26.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_D7_trial27.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_D8_trial28.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_E2_trial29.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_E3_trial30.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_E4_trial31.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_E5_trial32.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_E6_trial33.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_E7_trial34.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_E8_trial35.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_F2_trial36.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_F3_trial37.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_F4_trial38.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_F5_trial39.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_F6_trial40.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_F7_trial41.mat',
    'Rat Pharmacological Dataset/Rat_MEA2Pharma_DIV22_raster_3D_F8_trial42.mat',
}

for i_datafile = 1:length(datafiles)
datafilename = datafiles{i_datafile}

% read the data file and set the raster variable
raster = load(datafilename);
raster = raster.raster_3D;

%load parameters from a file

var_params = readmatrix("parameter_file.txt");

x_window = var_params(1);
y_window = var_params(2);
time_window = var_params(3);
n_noise_implementations = var_params(4);
num_epochs = var_params(5);
epoch_length = var_params(6);

max_time_lag = ceil(time_window/2);

size_snippet = 2*max_time_lag + epoch_length + 1;
total_length = size(raster,3);

n_epochs = floor(total_length / epoch_length );

if (num_epochs == 1)
    n_epochs = 1;
end

actual_contribution=cell(n_epochs,1);
actual_noise = cell(n_epochs,1);
an_ratio = cell(n_epochs,1);

init = 1;
post_end = 0;

for ii = init:n_epochs
    tic
    disp(ii)

    if (post_end + 1 + size_snippet) < size(raster,3)

        slice_start = post_end+1;

        pre_snippet_raster = raster(:,:,slice_start:slice_start+max_time_lag-1);
        pre_end = slice_start+max_time_lag-1;

        snippet_raster = raster(:,:, pre_end+1: pre_end+1+epoch_length);
        snip_end = pre_end+1+epoch_length;

        post_snippet_raster = raster(:,:, snip_end+1 : snip_end+1+max_time_lag-1);
        post_end = snip_end+1+max_time_lag-1;

        temp_snippet = cat(3,pre_snippet_raster, snippet_raster, post_snippet_raster);
        
        size(temp_snippet)

        tic
        [class_contribution]...
          = triple_correlation_class_contributions_3D_sp_wr(temp_snippet, x_window, ...
                        y_window, time_window);
        toc

        actual_contribution{ii} = class_contribution;

        toc
    end
end
str1 = 'COMPLETED_';
str2 = datafilename;

resultfile = append(str1,str2);
save(resultfile);
end