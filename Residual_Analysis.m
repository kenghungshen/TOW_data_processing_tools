%addpath('C:\Users\ks49639\Documents\MATLAB\ezc3d_matlab')
clear
clc
%% Load data
% Load FP data and store in a matrix with 6 columns (one for each channel) 

% Get the data from unload trial
    % Trial used to remove inertial force
Unload_file_name_all = ["Unload_SSOS.c3d" "Unload_FSOS.c3d"];

for i = 1:2
unload_file_name = convertStringsToChars(Unload_file_name_all(i));
Unload_OS  = ezc3dRead(unload_file_name); 

    %analog
    unload_analog_header = string(Unload_OS.parameters.ANALOG.LABELS.DATA);
        %find signal location handle in the matrix
    F1X =  find(unload_analog_header=="Force.Fx1");
    F1Y = find(unload_analog_header=="Force.Fx1");
    F1Z = find(unload_analog_header=="Force.Fx1");
    M1X = find(unload_analog_header=="Moment.Mx1");
    M1Y = find(unload_analog_header=="Moment.My1");
    M1Z = find(unload_analog_header=="Moment.Mz1");

    F2X =  find(unload_analog_header=="Force.Fx2");
    F2Y = find(unload_analog_header=="Force.Fx1");
    F2Z = find(unload_analog_header=="Force.Fx1");
    M2X = find(unload_analog_header=="Moment.Mx2");
    M2Y = find(unload_analog_header=="Moment.My2");
    M2Z = find(unload_analog_header=="Moment.Mz2");
    
    Unload_F1X_sig = Unload_OS.data.analogs(:,F1X);
    Unload_F1Y_sig =Unload_OS.data.analogs(:,F1Y);
    Unload_F1Z_sig = Unload_OS.data.analogs(:,F1Z);
    Unload_M1X_sig =Unload_OS.data.analogs(:,M1X);
    Unload_M1Y_sig = Unload_OS.data.analogs(:,M1Y);
    Unload_M1Z_sig =Unload_OS.data.analogs(:,M1Z);

    Unload_F2X_sig = Unload_OS.data.analogs(:,F2X);
    Unload_F2Y_sig =Unload_OS.data.analogs(:,F2Y);
    Unload_F2Z_sig = Unload_OS.data.analogs(:,F2Z);
    Unload_M2X_sig =Unload_OS.data.analogs(:,M2X);
    Unload_M2Y_sig = Unload_OS.data.analogs(:,M2Y);
    Unload_M2Z_sig =Unload_OS.data.analogs(:,M2Z);

    % sort as a matrix
    if contains(unload_file_name,'SS')
        Unload_SS_FP_sig = [Unload_F1X_sig, Unload_F1Y_sig, Unload_F1Z_sig, Unload_M1X_sig, Unload_M1Y_sig, Unload_M1Z_sig,...
            Unload_F2X_sig, Unload_F2Y_sig, Unload_F2Z_sig, Unload_M2X_sig, Unload_M2Y_sig, Unload_M2Z_sig];
    elseif contains(unload_file_name,'FS')
        Unload_FS_FP_sig = [Unload_F1X_sig, Unload_F1Y_sig, Unload_F1Z_sig, Unload_M1X_sig, Unload_M1Y_sig, Unload_M1Z_sig,...
            Unload_F2X_sig, Unload_F2Y_sig, Unload_F2Z_sig, Unload_M2X_sig, Unload_M2Y_sig, Unload_M2Z_sig];
    end
    %marker
    marker_header = string(Unload_OS.parameters.POINT.LABELS.DATA);
    LtTread1 = find(marker_header=="LeftTreadmill1");
    RtTread1 = find(marker_header=="RightTreadmill1");

    if contains(unload_file_name,'SS')
    Unload_SS_LtTread1X = squeeze(Unload_OS.data.points(1,LtTread1,:));
    Unload_SS_RtTread1X = squeeze(Unload_OS.data.points(1,RtTread1,:));
    elseif contains(unload_file_name,'FS')
    Unload_FS_LtTread1X = squeeze(Unload_OS.data.points(1,LtTread1,:));
    Unload_FS_RtTread1X = squeeze(Unload_OS.data.points(1,RtTread1,:));
    end

    clear Unload_OS
end

    Unload_FP_sig = {Unload_SS_FP_sig, Unload_FS_FP_sig};
    Unload_LtTread1X = {Unload_SS_LtTread1X, Unload_FS_LtTread1X};
    Unload_RtTread1X = {Unload_SS_RtTread1X, Unload_FS_RtTread1X};
    
%% Get the c3d file to be adjusted
trial_name = ["BaseSS01.c3d" "BaseSS02.c3d" "BaseFS01.c3d" "BaseFS02.c3d"...
     "BaseSSOS01.c3d" "BaseSSOS02.c3d" "BaseFSOS01.c3d" "BaseFSOS02.c3d"...
     "BaseOSSS01.c3d" "BaseOSFS01.c3d"];

 % type in the trial that we want to do residual analysis
 
trial = 6%1:length(trial_name) | 1- 4 = SS or FS, 5 - 8 = w/OS;
    file_name = convertStringsToChars(trial_name(trial))
    OS = contains(file_name,'OS0');
    % input into compensation function and sort as matrix for residual
    % analysis
    force_adjustment(file_name, OS,trial, Unload_LtTread1X, Unload_RtTread1X, Unload_FP_sig,Unload_file_name_all)

%% Residual Analysis
R_sig_col = [7:12];
L_sig_col = [1:6];

sig_col = R_sig_col;
%sig_col = L_sig_col;

for osci = 1:OS+1

% Call the input signal: "signal_data".
if osci == 1
signal_data = Original_FP(:,sig_col);
elseif osci ==2
signal_data = Compensated_FP(:,sig_col);
end

% Define a range of possible cutoff frequencies (in Hz)
cutoff_freqs = 1:1:100;

% for plotting and defining optimal regression range
ttl = ["FX", "FY", "FZ", "MX", "MY", "MZ"];
reg_range = [30:1:100];

% Set the sampling frequency (in Hz)
Fs = 2000;

% Define the filter order
filter_order = 4;

% Preallocate a matrix to store the residuals for each channel and cutoff frequency
residuals = zeros(length(signal_data), 6, length(cutoff_freqs));

% Loop over each channel and each cutoff frequency
for i = 1:6 % each channel
    for j = 1:length(cutoff_freqs)
        % lowpass with the current cutoff frequency and filter order
        [b, a] = butter(filter_order, cutoff_freqs(j)/(Fs/2), 'low');
        filtered_data = filtfilt(b, a, signal_data(:,i));
        
        % Compute the residuals by subtracting the filtered data from the original data
        residuals(:,i,j) = signal_data(:,i) - filtered_data;
    end
end

% Compute the RMS (root mean square) of the residuals for each channel and cutoff frequency
rms_res = sqrt(mean(residuals.^2, 1));

% Plot the RMS residuals for each channel and cutoff frequency

for i = 1:6
    figure(osci*2-1);
    sgtitle(file_name)
    subplot(2,3,i);
    plot(cutoff_freqs, squeeze(rms_res(1,i,:)));
    xlabel('Cutoff Frequency (Hz)');
    ylabel('RMS Residual');
    title(ttl(i));
    
    % fit regression line for the tail of residual line
    mdl = fitlm(reg_range,squeeze(rms_res(1,i,reg_range)));
    line([0,cutoff_freqs(end)],[mdl.Coefficients.Estimate(1),cutoff_freqs(end)*mdl.Coefficients.Estimate(2)+mdl.Coefficients.Estimate(1)],...
        'Color','r','LineStyle','--')
    % plot interept and find intersection with residual line
    yline(mdl.Coefficients.Estimate(1),'Color','g','LineStyle','--')
    [~,idx] = min(abs(squeeze(rms_res(1,i,:)) - mdl.Coefficients.Estimate(1)));
    xline(idx,'Color','k','LineStyle','--')
    text(mean(reg_range),mean(squeeze(rms_res(1,i,:))),['optimal cutoff = ' num2str(idx)]);
    
    figure(osci*2);
    sgtitle('Filtered with optimal cutoff')
    %idx = 15; % assign specific filter cutoff freq that's different from
    %the finding of residual analysis
    subplot(2,3,i);
    plot(signal_data(:,i) - squeeze(residuals(:,i,idx)));
    title([ttl(i) 'filtered at' num2str(idx) 'Hz']);
    
end
end

%% Compensation algorithm
function force_adjustment(file_name,OS,trial, Unload_LtTread1X, Unload_RtTread1X, Unload_FP_sig,Unload_file_name_all)

if isfile(file_name)

 disp(['Processing file: ',file_name])

 if OS == 1
     disp('Force Compensation = True ')
 elseif OS == 0
     disp('Force Compensation = False ')
 end

TOW  = ezc3dRead(file_name); 
disp(['analog sampling rate = ' num2str(TOW.header.analogs.frameRate)]);

% Get the force data to be adjusted from file TOW
    %analog
TOW_analog_header = string(TOW.parameters.ANALOG.LABELS.DATA);
F1X =  find(TOW_analog_header=="Force.Fx1");
F1Y = find(TOW_analog_header=="Force.Fx1");
F1Z = find(TOW_analog_header=="Force.Fx1");
M1X = find(TOW_analog_header=="Moment.Mx1");
M1Y = find(TOW_analog_header=="Moment.My1");
M1Z = find(TOW_analog_header=="Moment.Mz1");

F2X =  find(TOW_analog_header=="Force.Fx2");
F2Y = find(TOW_analog_header=="Force.Fx1");
F2Z = find(TOW_analog_header=="Force.Fx1");
M2X = find(TOW_analog_header=="Moment.Mx2");
M2Y = find(TOW_analog_header=="Moment.My2");
M2Z = find(TOW_analog_header=="Moment.Mz2");


TOW_F1X_sig = TOW.data.analogs(:,F1X);
TOW_F1Y_sig =TOW.data.analogs(:,F1Y);
TOW_F1Z_sig = TOW.data.analogs(:,F1Z);
TOW_M1X_sig =TOW.data.analogs(:,M1X);
TOW_M1Y_sig = TOW.data.analogs(:,M1Y);
TOW_M1Z_sig =TOW.data.analogs(:,M1Z);

TOW_F2X_sig = TOW.data.analogs(:,F2X);
TOW_F2Y_sig =TOW.data.analogs(:,F2Y);
TOW_F2Z_sig = TOW.data.analogs(:,F2Z);
TOW_M2X_sig =TOW.data.analogs(:,M2X);
TOW_M2Y_sig = TOW.data.analogs(:,M2Y);
TOW_M2Z_sig =TOW.data.analogs(:,M2Z);

TOW_FP_sig = [TOW_F1X_sig, TOW_F1Y_sig, TOW_F1Z_sig, TOW_M1X_sig, TOW_M1Y_sig, TOW_M1Z_sig,...
    TOW_F2X_sig, TOW_F2Y_sig, TOW_F2Z_sig, TOW_M2X_sig, TOW_M2Y_sig, TOW_M2Z_sig];

if OS == 1 % go into compensation algorithm
    if contains(file_name,'SS')
        speed_handle = 1;
    elseif contains(file_name,'FS')
        speed_handle = 2;
    end
    
    disp(['Compensating using file: ', convertStringsToChars(Unload_file_name_all(speed_handle))])
    
    %marker
TOW_marker_header = string(TOW.parameters.POINT.LABELS.DATA);
LtTread1 = find(TOW_marker_header=="LeftTreadmill1");
RtTread1 = find(TOW_marker_header=="RightTreadmill1");

TOW_LtTread1X = squeeze(TOW.data.points(1,LtTread1,:));
TOW_RtTread1X = squeeze(TOW.data.points(1,RtTread1,:));

% Force Compensate
    
% sync using treadmill marker
    % Find the first to cross zero after demean
    % Using right treadmill only for now
    % speed handle determine which unload trial to use
    unload_tread_handle = Unload_RtTread1X{1,speed_handle};
    unload_start = find((unload_tread_handle(1:end-1)-mean(unload_tread_handle)).*(unload_tread_handle(2:end)-mean(unload_tread_handle))<0,4);
    TOW_start = find((TOW_RtTread1X(1:end-1)-mean(TOW_RtTread1X)).*(TOW_RtTread1X(2:end)-mean(TOW_RtTread1X))<0,1);
    unload_start = unload_start(3:4); % remove first 2 in unload trial - make sure the instant happens in TOW earlier then unload
    
    if diff(TOW_RtTread1X(TOW_start:TOW_start+1))*diff(unload_tread_handle(unload_start(1):unload_start(1)+1))>0
    unload_start = unload_start(1);
    else
    unload_start = unload_start(2);
    end

    
    % Truncate unload trial to Sync Kinematics
    unload_truncate_kinematic = [unload_start - TOW_start + 1 : unload_start - TOW_start + 1 + length(TOW_RtTread1X)-1];
        
        % Transform kinematic frames to kinetic
    unload_truncate_kinetic_start = unload_truncate_kinematic(1)*20-19;
    unload_truncate_kinetic = [unload_truncate_kinetic_start : unload_truncate_kinetic_start+length(TOW_FP_sig)-1]';
    
    % Compensate
    TOW_FP_compensated = TOW_FP_sig - Unload_FP_sig{1,speed_handle}(unload_truncate_kinetic);
    disp('Compensation Completed')
    

    % Original signal: TOW_FP_sig
    % compensate signal: TOW_FP_compensated
    
% Quality Check
    % Check the sync of marker
    % L
    
    for i = 1:6
        if i ==1
            fig_title = " Fx";
        elseif i==2
            fig_title = " Fy";
        elseif i==3
            fig_title = " Fz";
        elseif i==4
            fig_title = " Mx";
        elseif i==5
            fig_title = " My";
        elseif i==6
            fig_title = " Mz";
        end
        
        figure(trial*10 + i)
        
    subplot(2,4,1)
    plot(TOW_LtTread1X(1:4000)-mean(TOW_LtTread1X),'b'); hold on
    plot(Unload_LtTread1X{1,speed_handle}(unload_truncate_kinematic)-mean(Unload_LtTread1X{1,speed_handle}),'r');
    title('Matching the marker trajectory (Lt Treadmill)')
    
    subplot(2,4,5)
    plot(TOW_FP_sig(:,i),'b'); hold on
    plot(Unload_FP_sig{1,speed_handle}(unload_truncate_kinetic,i),'r'); 
    title('Matching the force')
    legend(strcat('TOW', fig_title),strcat('Unloaded', fig_title))
    
    plot_min = min([TOW_FP_sig(:,i) ; Unload_FP_sig{1,speed_handle}(unload_truncate_kinetic,i)]); 
    plot_max = max([TOW_FP_sig(:,i) ; Unload_FP_sig{1,speed_handle}(unload_truncate_kinetic,i)]); 
    ylim([plot_min,plot_max])
    
    subplot(2,4,2)
    plot(TOW_FP_sig(:,i),'b'); hold on
    title(strcat('original',fig_title))
    ylim([plot_min,plot_max])
    
    subplot(2,4,6)
    plot(TOW_FP_compensated(:,i),'b');
    title(strcat('compensated', fig_title))
    ylim([plot_min,plot_max])
    
    
    % R
    ii = i+6;
    subplot(2,4,3)
    plot(TOW_RtTread1X(1:4000)-mean(TOW_RtTread1X),'b'); hold on
    plot(unload_tread_handle(unload_truncate_kinematic)-mean(unload_tread_handle),'r');
    title('Matching the marker trajectory (Rt Treadmill)')
    
    subplot(2,4,7)
    plot(TOW_FP_sig(:,ii),'b'); hold on
    plot(Unload_FP_sig{1,speed_handle}(unload_truncate_kinetic,ii),'r'); 
    title('Matching the force')
    legend(strcat('TOW', fig_title),strcat('Unloaded', fig_title))
    
    plot_min = min([TOW_FP_sig(:,ii) ; Unload_FP_sig{1,speed_handle}(unload_truncate_kinetic,ii)]); 
    plot_max = max([TOW_FP_sig(:,ii) ; Unload_FP_sig{1,speed_handle}(unload_truncate_kinetic,ii)]); 
    ylim([plot_min,plot_max])

    
    subplot(2,4,4)
    plot(TOW_FP_sig(:,ii),'b'); hold on
    title(strcat('original',fig_title))
    ylim([plot_min,plot_max])
    
    subplot(2,4,8)
    plot(TOW_FP_compensated(:,ii),'b');
    title(strcat('compensated', fig_title))
    ylim([plot_min,plot_max])

    end
    
end
end

% load into global workspace
assignin('base', 'Original_FP', TOW_FP_sig);
if OS ==1
assignin('base', 'Compensated_FP', TOW_FP_compensated);
end

end