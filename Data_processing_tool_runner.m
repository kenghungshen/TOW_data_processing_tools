addpath('C:\Users\ks49639\Documents\MATLAB\ezc3d_matlab')
clear
clc

% recording the command window output
diary off
filename = 'Matlab_output.txt';
diary(filename)

%% Filter parameter
force_filter_cutoff = 15

EMG_bandpass_cutoff = [16 500]
EMG_lowpass_cutoff = 20

%% Inertial Compensation channels
%compensate_channel = "Fx"
compensate_channel = "All"

Comp_QC_plot = "All";
%Comp_QC_plot = "compensated"
%% Trial names
% Static
static_file_name = 'Static.c3d'

weight_scale_data = 100;
weight_scale_session = ["NA"]

% vGRFCadenceTest
% vCT_steps = 20;
% vCT_file_name = ["vGRFCadenceTest_SS.c3d" "vGRFCadenceTest_FS.c3d"]

% Force processing
Unload_file_name_all = ["Unload_SSOS.c3d" "Unload_FSOS.c3d"]

trial_name = ["BaseSS01.c3d" "BaseSS02.c3d" "BaseFS01.c3d" "BaseFS02.c3d"...
     "BaseSSOS01.c3d" "BaseSSOS02.c3d" "BaseFSOS01.c3d" "BaseFSOS02.c3d"...
     "BaseOSSS01.c3d" "BaseOSFS01.c3d"]%...
     %"PostSSnew01.c3d" "PostSSnew02.c3d" "PostFSnew01.c3d" "PostFSnew02.c3d"]
 
% EMG processing
Fadjusted_trial_name = strcat("Fadjusted_", trial_name);

%% EMG sensor number
LVL_sensor = 'Sensor 1.EMG1'
LGM_sensor = 'Sensor 3.EMG3'
LTFL_sensor = 'Sensor 4.EMG4'

RVL_sensor = 'Sensor 5.EMG5'
RGM_sensor = 'Sensor 6.EMG6'
RTFL_sensor = 'Sensor 8.EMG8'
%% Running scripts
Inertial_Compensation
c3d_EMG_filtering
%vGRFCadenceTest_validation
%vCT: 4/12/2023 found figure off, must be bug in the code. Will fix later. Skip for now and manually check during data collection.
FP_LAB_Match
Static_BW_validation
% Static_BW_validation: 8/19/2023 some post-training data has no weght scale data, added this module to validate using static trail

%%
diary off