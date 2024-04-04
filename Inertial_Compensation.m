%addpath('C:\Users\ks49639\Documents\MATLAB\ezc3d_matlab')

%%
%clc
%clear

%% Force filtering parameters
% force_filter_cutoff = 15;

%% Compensate channels
% compensate_channel = "Fx"
%compensate_channel = "All"
%% Get the data from unload trial
% Trial used to remove inertial force
% 
% Unload_file_name_all = ["Unload_SSOS.c3d" "Unload_FSOS.c3d"];

for i = 1:2
unload_file_name = convertStringsToChars(Unload_file_name_all(i));
Unload_OS  = ezc3dRead(unload_file_name); 

    %analog
    unload_analog_header = string(Unload_OS.parameters.ANALOG.LABELS.DATA);
        %find signal location handle in the matrix
    F1X =  find(unload_analog_header=="Force.Fx1");
    F1Y = find(unload_analog_header=="Force.Fy1");
    F1Z = find(unload_analog_header=="Force.Fz1");
    M1X = find(unload_analog_header=="Moment.Mx1");
    M1Y = find(unload_analog_header=="Moment.My1");
    M1Z = find(unload_analog_header=="Moment.Mz1");

    F2X =  find(unload_analog_header=="Force.Fx2");
    F2Y = find(unload_analog_header=="Force.Fy2");
    F2Z = find(unload_analog_header=="Force.Fz2");
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
    
    disp('Finish Processing the Unload trial, now start working on the walking trials.')
    disp(' ')
    
    
%% Get the c3d file to be adjusted
% trial_name = ["BaseSS01.c3d" "BaseSS02.c3d" "BaseFS01.c3d" "BaseFS02.c3d"...
%      "BaseSSOS01.c3d" "BaseSSOS02.c3d" "BaseFSOS01.c3d" "BaseFSOS02.c3d"...
%      "BaseOSSS01.c3d" "BaseOSFS01.c3d"];

for trial = 1:length(trial_name)
    file_name = convertStringsToChars(trial_name(trial));
    OS = contains(file_name,'OS0');
    force_adjustment(file_name, OS, force_filter_cutoff, Unload_LtTread1X, Unload_RtTread1X, Unload_FP_sig,Unload_file_name_all,compensate_channel, Comp_QC_plot)
end
 
%% Function to adjust force signal (compensate if needed, then low pass filter)
function force_adjustment(file_name,OS, force_filter_cutoff, Unload_LtTread1X, Unload_RtTread1X, Unload_FP_sig,Unload_file_name_all,compensate_channel, Comp_QC_plot)

if isfile(file_name)
    
 disp(['Processing file: ',file_name])

 if OS == 1
     disp('Force Compensation = True ')
 elseif OS == 0
     disp('Force Compensation = False ')
 end

TOW  = ezc3dRead(file_name); 

% The force data to be adjusted
% TOW
    %analog
TOW_analog_header = string(TOW.parameters.ANALOG.LABELS.DATA);
F1X =  find(TOW_analog_header=="Force.Fx1");
F1Y = find(TOW_analog_header=="Force.Fy1");
F1Z = find(TOW_analog_header=="Force.Fz1");
M1X = find(TOW_analog_header=="Moment.Mx1");
M1Y = find(TOW_analog_header=="Moment.My1");
M1Z = find(TOW_analog_header=="Moment.Mz1");

F2X =  find(TOW_analog_header=="Force.Fx2");
F2Y = find(TOW_analog_header=="Force.Fy2");
F2Z = find(TOW_analog_header=="Force.Fz2");
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

if OS == 1
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
    %% Force Compensate
    
% sync using treadmill marker
    % Find the first to cross zero after demean
    % Using "right" treadmill to sync; plot both R and L to check
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
    ana_Fs = TOW.header.analogs.frameRate;
    marker_Fs = TOW.header.points.frameRate;
    Fs_kine_kinema = ana_Fs/marker_Fs;
    disp(['analog Fs = ' num2str(ana_Fs)])
    disp(['Kinematic Fs = ' num2str(marker_Fs)])
    disp(['Fs analog/kinematics = ' num2str(Fs_kine_kinema)])
    
    unload_truncate_kinetic_start = unload_truncate_kinematic(1)*Fs_kine_kinema-(Fs_kine_kinema-1);
    unload_truncate_kinetic = [unload_truncate_kinetic_start : unload_truncate_kinetic_start+length(TOW_FP_sig)-1]';
    
    % Compensate
    
    if compensate_channel =="Fx"
        Unload_FP_sig{1,speed_handle}(unload_truncate_kinetic,[2:6, 8:12]) = 0;
        if Comp_QC_plot == "All"
        plot_channel = 1:6;
        else
            plot_channel = 1;
        end
        disp('Only compensate Fx, Removed Unload Force channels in Unload FP Sig')
    elseif compensate_channel == "All"
        plot_channel = 1:6;
        disp('Compensate all channels. Unload FP Sig not modified')
    end

    TOW_FP_compensated = TOW_FP_sig - Unload_FP_sig{1,speed_handle}(unload_truncate_kinetic,:);
    disp('Compensation Completed')
    

    % Filter force with and without compensation for quality check
    % %low pass filter
     f_analog = TOW.header.analogs.frameRate; % sampling rate
     f_analog = f_analog/2; % Nyquest frequency
     f_cutoff = force_filter_cutoff;  % lowpass threshold
     disp('Now lowpass filtering the FP signal')
     disp(['Analog sampling rate/2 = ', num2str(f_analog)]);
     disp(['Cutoff frequency= ', num2str(f_cutoff)]);
     
     order = 4; % 4th order
     [c,e] = butter(order,f_cutoff/f_analog,'low');
 
     TOW_FP_original_filt=filtfilt(c,e,TOW_FP_sig); % filter lowpass 30Hz
     TOW_FP_compensated_filt=filtfilt(c,e,TOW_FP_compensated); % filter lowpass 30Hz

    
    %% Quality Check
    % Check the sync of marker
    
    for i = plot_channel%now only plot Fx. If want to plot all 6 channels: change to i=1:6
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
        
    figure(i);
    subplot(4,3,1)
    plot(TOW_LtTread1X(1:4000)-mean(TOW_LtTread1X),'b'); hold on
    plot(Unload_LtTread1X{1,speed_handle}(unload_truncate_kinematic)-mean(Unload_LtTread1X{1,speed_handle}),'r');
    title('Matching the marker trajectory (Lt Treadmill)')
    
    subplot(4,3,4)
    plot(TOW_FP_sig(:,i),'b'); hold on
    plot(Unload_FP_sig{1,speed_handle}(unload_truncate_kinetic,i),'r'); 
    title('Matching the force')
    legend(strcat('TOW', fig_title),strcat('Unloaded', fig_title))
    
    plot_min = min([TOW_FP_sig(:,i) ; Unload_FP_sig{1,speed_handle}(unload_truncate_kinetic,i)]); 
    plot_max = max([TOW_FP_sig(:,i) ; Unload_FP_sig{1,speed_handle}(unload_truncate_kinetic,i)]); 
    ylim([plot_min,plot_max])
    
    subplot(4,3,2)
    plot(TOW_FP_sig(:,i),'b'); hold on
    title(strcat('original',fig_title))
    ylim([plot_min,plot_max])
    
    subplot(4,3,5)
    plot(TOW_FP_compensated(:,i),'b');
    title(strcat('compensated', fig_title))
    ylim([plot_min,plot_max])
    
    subplot(4,3,3)
    plot(TOW_FP_original_filt(:,i),'b'); hold on
    title(strcat('original filtered',fig_title))
    ylim([plot_min,plot_max])
    
    subplot(4,3,6)
    plot(TOW_FP_compensated_filt(:,i),'b');
    title(strcat('compensated filtered',fig_title))
    ylim([plot_min,plot_max])
    
    % R
    ii = i+6;
    subplot(4,3,7)
    plot(TOW_RtTread1X(1:4000)-mean(TOW_RtTread1X),'b'); hold on
    plot(unload_tread_handle(unload_truncate_kinematic)-mean(unload_tread_handle),'r');
    title('Matching the marker trajectory (Rt Treadmill)')
    
    subplot(4,3,10)
    plot(TOW_FP_sig(:,ii),'b'); hold on
    plot(Unload_FP_sig{1,speed_handle}(unload_truncate_kinetic,ii),'r'); 
    title('Matching the force')
    legend(strcat('TOW', fig_title),strcat('Unloaded', fig_title))
    
    plot_min = min([TOW_FP_sig(:,ii) ; Unload_FP_sig{1,speed_handle}(unload_truncate_kinetic,ii)]); 
    plot_max = max([TOW_FP_sig(:,ii) ; Unload_FP_sig{1,speed_handle}(unload_truncate_kinetic,ii)]); 
    ylim([plot_min,plot_max])

    
    subplot(4,3,8)
    plot(TOW_FP_sig(:,ii),'b'); hold on
    title(strcat('original',fig_title))
    ylim([plot_min,plot_max])
    
    subplot(4,3,11)
    plot(TOW_FP_compensated(:,ii),'b');
    title(strcat('compensated', fig_title))
    ylim([plot_min,plot_max])
    
    subplot(4,3,9)
    plot(TOW_FP_original_filt(:,ii),'b'); hold on
    title(strcat('original filtered',fig_title))
    ylim([plot_min,plot_max])
    
    subplot(4,3,12)
    plot(TOW_FP_compensated_filt(:,ii),'b');
    title(strcat('compensated filtered',fig_title))
    ylim([plot_min,plot_max])
    end
    
    uiwait(gcf);
    
elseif OS == 0

     f_analog = TOW.header.analogs.frameRate; % sampling rate
     f_analog = f_analog/2; % Nyquest frequency
     f_cutoff = force_filter_cutoff;  % lowpass threshold
     
     disp('Now lowpass filtering the FP signal')
     disp(['Analog sampling rate/2 = ', num2str(f_analog)]);
     disp(['Cutoff frequency= ', num2str(f_cutoff)]);
     
     order = 4; % 4th order
     [c,e] = butter(order,f_cutoff/f_analog,'low');
 
     TOW_FP_filt=filtfilt(c,e,TOW_FP_sig); % filter lowpass 30Hz
     
end
    %% Write c3d

if OS == 1
    TOW.data.analogs(:,F1X)=TOW_FP_compensated_filt(:,1);
    TOW.data.analogs(:,F1Y)=TOW_FP_compensated_filt(:,2);
    TOW.data.analogs(:,F1Z)=TOW_FP_compensated_filt(:,3);
    TOW.data.analogs(:,M1X)=TOW_FP_compensated_filt(:,4);
    TOW.data.analogs(:,M1Y)=TOW_FP_compensated_filt(:,5);
    TOW.data.analogs(:,M1Z)=TOW_FP_compensated_filt(:,6);

    TOW.data.analogs(:,F2X)=TOW_FP_compensated_filt(:,7);
    TOW.data.analogs(:,F2Y)=TOW_FP_compensated_filt(:,8);
    TOW.data.analogs(:,F2Z)=TOW_FP_compensated_filt(:,9);
    TOW.data.analogs(:,M2X)=TOW_FP_compensated_filt(:,10);
    TOW.data.analogs(:,M2Y)=TOW_FP_compensated_filt(:,11);
    TOW.data.analogs(:,M2Z)=TOW_FP_compensated_filt(:,12);
    
elseif OS == 0
    TOW.data.analogs(:,F1X)=TOW_FP_filt(:,1);
    TOW.data.analogs(:,F1Y)=TOW_FP_filt(:,2);
    TOW.data.analogs(:,F1Z)=TOW_FP_filt(:,3);
    TOW.data.analogs(:,M1X)=TOW_FP_filt(:,4);
    TOW.data.analogs(:,M1Y)=TOW_FP_filt(:,5);
    TOW.data.analogs(:,M1Z)=TOW_FP_filt(:,6);

    TOW.data.analogs(:,F2X)=TOW_FP_filt(:,7);
    TOW.data.analogs(:,F2Y)=TOW_FP_filt(:,8);
    TOW.data.analogs(:,F2Z)=TOW_FP_filt(:,9);
    TOW.data.analogs(:,M2X)=TOW_FP_filt(:,10);
    TOW.data.analogs(:,M2Y)=TOW_FP_filt(:,11);
    TOW.data.analogs(:,M2Z)=TOW_FP_filt(:,12);
end

     new_file_name = ['Fadjusted_',file_name];
     disp(['Writing File: ',new_file_name]);
     ezc3dWrite(new_file_name, TOW);
     disp('Done ------------------------------------------------------------------------------------------------------');
     disp(' ')
     
end
end

%%

% % sync using time
% t_marker = [0:1/ICT.header.points.frameRate:ICT.header.points.lastFrame/ICT.header.points.frameRate-1/ICT.header.points.frameRate];
% t_analog = [0:1/ICT.header.analogs.frameRate:ICT.header.analogs.lastFrame/ICT.header.analogs.frameRate-1/ICT.header.analogs.frameRate];