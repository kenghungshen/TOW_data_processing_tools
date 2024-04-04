%addpath('C:\Users\ks49639\Documents\MATLAB\ezc3d_matlab')

%clc
%clear
%% EMG filtering parameters
% bandpass - de-mean - rectify - lowpass

% EMG_bandpass_cutoff = [16 500];
% EMG_lowpass_cutoff = 20;

%% Link muscle to sensor number / analog name in c3d
% LVL_sensor = 'Sensor 1.EMG1';
% LGM_sensor = 'Sensor 3.EMG3';
% LTFL_sensor = 'Sensor 4.EMG4';
% 
% RVL_sensor = 'Sensor 5.EMG5';
% RGM_sensor = 'Sensor 6.EMG6';
% RTFL_sensor = 'Sensor 7.EMG7';

%% Get the data from MVIC trials
% get data from each muscle
for side = ["R", "L"]
    for joint=["Hip", "Knee"]
        for trial = 1:2
            file_name = ['MVIC-',convertStringsToChars(side),convertStringsToChars(joint),'0',num2str(trial),'.c3d'];
            disp(['Now processing' file_name]);
            current_file = ezc3dRead(file_name);
            
            f_analog = current_file.header.analogs.frameRate;
            analog_header = string(current_file.parameters.ANALOG.LABELS.DATA);
            
            if  joint == "Knee"
                if side == "R"
                    loc_RVL =  find(analog_header==RVL_sensor);          
                    MVIC_RVL_temp = current_file.data.analogs(:,loc_RVL);
                    filtered_MVIC_RVL_temp = EMG_filter(f_analog,EMG_bandpass_cutoff,EMG_lowpass_cutoff,MVIC_RVL_temp);
                    MVIC_RVL_temp = abs(MVIC_RVL_temp - mean(MVIC_RVL_temp)); % de-meaned and rectified raw signal for plotting
                    if trial == 1
                        MVIC_RVL = MVIC_RVL_temp; % raw
                        filtered_MVIC_RVL =  filtered_MVIC_RVL_temp; % filtered
                    elseif trial == 2
                        MVIC_RVL = [MVIC_RVL; MVIC_RVL_temp];
                        filtered_MVIC_RVL =  [filtered_MVIC_RVL; filtered_MVIC_RVL_temp];
                    end
                                        
                elseif side == "L"
                    loc_LVL =  find(analog_header==LVL_sensor);          
                    MVIC_LVL_temp = current_file.data.analogs(:,loc_LVL);
                    filtered_MVIC_LVL_temp = EMG_filter(f_analog,EMG_bandpass_cutoff,EMG_lowpass_cutoff,MVIC_LVL_temp);
                    MVIC_LVL_temp = abs(MVIC_LVL_temp - mean(MVIC_LVL_temp));
                    if trial == 1
                        MVIC_LVL = MVIC_LVL_temp;
                        filtered_MVIC_LVL =  filtered_MVIC_LVL_temp;
                    elseif trial == 2
                        MVIC_LVL = [MVIC_LVL; MVIC_LVL_temp];
                        filtered_MVIC_LVL =  [filtered_MVIC_LVL; filtered_MVIC_LVL_temp];
                    end
                end
                
                    
            elseif joint == "Hip"
                if side == "R"
                    loc_RGM =  find(analog_header==RGM_sensor);          
                    MVIC_RGM_temp = current_file.data.analogs(:,loc_RGM);
                    filtered_MVIC_RGM_temp = EMG_filter(f_analog,EMG_bandpass_cutoff,EMG_lowpass_cutoff,MVIC_RGM_temp);
                    MVIC_RGM_temp = abs(MVIC_RGM_temp - mean(MVIC_RGM_temp));
                    
                    loc_RTFL =  find(analog_header==RTFL_sensor);          
                    MVIC_RTFL_temp = current_file.data.analogs(:,loc_RTFL);
                    filtered_MVIC_RTFL_temp = EMG_filter(f_analog,EMG_bandpass_cutoff,EMG_lowpass_cutoff,MVIC_RTFL_temp);
                    MVIC_RTFL_temp = abs(MVIC_RTFL_temp - mean(MVIC_RTFL_temp));
                    
                    if trial == 1
                        MVIC_RGM = MVIC_RGM_temp;
                        filtered_MVIC_RGM =  filtered_MVIC_RGM_temp;
                        MVIC_RTFL = MVIC_RTFL_temp;
                        filtered_MVIC_RTFL =  filtered_MVIC_RTFL_temp;
                    elseif trial == 2
                        MVIC_RGM = [MVIC_RGM; MVIC_RGM_temp];
                        filtered_MVIC_RGM =  [filtered_MVIC_RGM; filtered_MVIC_RGM_temp];
                        MVIC_RTFL = [MVIC_RTFL; MVIC_RTFL_temp];
                        filtered_MVIC_RTFL =  [filtered_MVIC_RTFL; filtered_MVIC_RTFL_temp];
                    end
                    
                elseif side == "L"
                    loc_LGM =  find(analog_header==LGM_sensor);          
                    MVIC_LGM_temp = current_file.data.analogs(:,loc_LGM);
                    filtered_MVIC_LGM_temp = EMG_filter(f_analog,EMG_bandpass_cutoff,EMG_lowpass_cutoff,MVIC_LGM_temp);
                    MVIC_LGM_temp = abs(MVIC_LGM_temp - mean(MVIC_LGM_temp));
                    
                    loc_LTFL =  find(analog_header==LTFL_sensor);          
                    MVIC_LTFL_temp = current_file.data.analogs(:,loc_LTFL);
                    filtered_MVIC_LTFL_temp = EMG_filter(f_analog,EMG_bandpass_cutoff,EMG_lowpass_cutoff,MVIC_LTFL_temp);        
                    MVIC_LTFL_temp = abs(MVIC_LTFL_temp - mean(MVIC_LTFL_temp));
                    
                    if trial == 1
                        MVIC_LGM = MVIC_LGM_temp;
                        filtered_MVIC_LGM =  filtered_MVIC_LGM_temp;
                        MVIC_LTFL = MVIC_LTFL_temp;
                        filtered_MVIC_LTFL =  filtered_MVIC_LTFL_temp;
                    elseif trial == 2
                        MVIC_LGM = [MVIC_LGM; MVIC_LGM_temp];
                        filtered_MVIC_LGM =  [filtered_MVIC_LGM; filtered_MVIC_LGM_temp];
                        MVIC_LTFL = [MVIC_LTFL; MVIC_LTFL_temp];
                        filtered_MVIC_LTFL =  [filtered_MVIC_LTFL; filtered_MVIC_LTFL_temp];
                    end                   
                end                
            end
            disp('Done --------------------------------------');
        end
    end
end

% make sure the order are the same for all 3 matrices
MVIC_col_name = ["LVL", "LGM", "LTFL", "RVL", "RGM", "RTFL"];
MVIC = {MVIC_LVL, MVIC_LGM, MVIC_LTFL, MVIC_RVL, MVIC_RGM, MVIC_RTFL};
filtered_MVIC = {filtered_MVIC_LVL, filtered_MVIC_LGM, filtered_MVIC_LTFL, filtered_MVIC_RVL, filtered_MVIC_RGM, filtered_MVIC_RTFL};

fig(1) = figure('Name','Filter validate');
for i = 1:6
    subplot(2,3,i)
    plot(MVIC{i},'Color',[0.5, 0.5, 0.5, 0.2]); hold on
    plot(filtered_MVIC{i},'b');
    title(MVIC_col_name(i))
end
%% Feed into GUI for manual confirm/selection

MVIC_GUI(filtered_MVIC,fig, MVIC_col_name)
uiwait(gcf);

%% Get data from gait trials - filter and normalize to MVIC

% Fadjusted_trial_name = ["Fadjusted_BaseSS01.c3d" "Fadjusted_BaseSS02.c3d" "Fadjusted_BaseFS01.c3d" "Fadjusted_BaseFS02.c3d"...
%      "Fadjusted_BaseSSOS01.c3d" "Fadjusted_BaseSSOS02.c3d" "Fadjusted_BaseFSOS01.c3d" "Fadjusted_BaseFSOS02.c3d"...
%      "Fadjusted_BaseOSSS01.c3d" "Fadjusted_BaseOSFS01.c3d"];

for trial = 1:length(Fadjusted_trial_name)
            file_name = convertStringsToChars(Fadjusted_trial_name(trial));
            c3d_EMG_filter(file_name, LVL_sensor, LGM_sensor, LTFL_sensor, RVL_sensor, RGM_sensor, RTFL_sensor, max_MVIC,EMG_bandpass_cutoff,EMG_lowpass_cutoff)
end

%% function for EMG process in gait trials: filter + normalization and write new .c3d
function c3d_EMG_filter(file_name, LVL_sensor, LGM_sensor, LTFL_sensor, RVL_sensor, RGM_sensor, RTFL_sensor, max_MVIC, EMG_bandpass_cutoff,EMG_lowpass_cutoff)

if isfile(file_name)
disp(['Now procressing file: ', file_name]);
            
            current_file = ezc3dRead(file_name);
            
            f_analog = current_file.header.analogs.frameRate;
            analog_header = string(current_file.parameters.ANALOG.LABELS.DATA);
            
            loc_LVL = find(analog_header==LVL_sensor);          
            loc_LGM = find(analog_header==LGM_sensor);  
            loc_LTFL = find(analog_header==LTFL_sensor);  
            loc_RVL = find(analog_header==RVL_sensor);  
            loc_RGM = find(analog_header==RGM_sensor);  
            loc_RTFL = find(analog_header==RTFL_sensor);
            
            LVL = current_file.data.analogs(:,loc_LVL);
            LGM = current_file.data.analogs(:,loc_LGM);
            LTFL = current_file.data.analogs(:,loc_LTFL);
            RVL = current_file.data.analogs(:,loc_RVL);
            RGM = current_file.data.analogs(:,loc_RGM);
            RTFL = current_file.data.analogs(:,loc_RTFL);
            
            emg = [LVL, LGM, LTFL, RVL, RGM, RTFL];
            filtered_emg = EMG_filter(f_analog,EMG_bandpass_cutoff,EMG_lowpass_cutoff,emg);
            filtered_emg_norm = filtered_emg./max_MVIC;
            
            % re-write
            current_file.parameters.ANALOG.LABELS.DATA(loc_LVL)={'LVL'};
            current_file.parameters.ANALOG.LABELS.DATA(loc_LGM)={'LGM'};
            current_file.parameters.ANALOG.LABELS.DATA(loc_LTFL)={'LTFL'};
            current_file.parameters.ANALOG.LABELS.DATA(loc_RVL)={'RVL'};
            current_file.parameters.ANALOG.LABELS.DATA(loc_RGM)={'RGM'};
            current_file.parameters.ANALOG.LABELS.DATA(loc_RTFL)={'RTFL'};
            
            current_file.data.analogs(:,loc_LVL) = filtered_emg_norm(:,1);
            current_file.data.analogs(:,loc_LGM) = filtered_emg_norm(:,2);
            current_file.data.analogs(:,loc_LTFL) = filtered_emg_norm(:,3);
            current_file.data.analogs(:,loc_RVL) = filtered_emg_norm(:,4);
            current_file.data.analogs(:,loc_RGM) = filtered_emg_norm(:,5);
            current_file.data.analogs(:,loc_RTFL) = filtered_emg_norm(:,6);
            
            disp('EMG signal processed');
            new_file_name = ['EMGfiltered_',file_name];
            disp(['Writing File: ',new_file_name]);
            ezc3dWrite(new_file_name, current_file);
            disp('Done ---------------------------------------------------------------------------------------------------------------');
            
end
end

    %% function for EMG filtering
    function filtered_emg = EMG_filter(f_analog,EMG_bandpass_cutoff,EMG_lowpass_cutoff,emg)
     
     f_analog = f_analog/2; % Nyquest frequency
     
     f_cutoff = EMG_bandpass_cutoff;  % bandpass threshold
     order = 4; % 4th order
     [c,e] = butter(order,f_cutoff/f_analog,'bandpass');
    
     disp(['Analog sampling rate/2 = ', num2str(f_analog)]);
     disp(['Bandpass filtering the EMG signal. Cutoff frequency= ', num2str(f_cutoff(1)),'~',num2str(f_cutoff(2))]);
     emg_h=filtfilt(c,e,emg); % filter banpass 16-500
     
     emg_nh=abs(emg_h-mean(emg_h)); % de-mean + rectify
     disp('EMG de-meaned and rectified')
     
      f_cutoff = EMG_lowpass_cutoff;
      order = 4;
      [f,g] = butter(order,f_cutoff/f_analog,'low'); 
      
      disp(['Lowpass filtering the EMG signal. Cutoff frequency= ', num2str(f_cutoff)]);
      filtered_emg=filtfilt(f,g,emg_nh); % low pass 50 after de-mean + rectify           
           
    end
 
%% function for interactive MVIC selection and storage
function MVIC_GUI(filtered_MVIC,fig, MVIC_col_name)
% Create figure and axes for subplots
fig(2) = figure('Name','MVIC EMG');

for i = 1:6
    ax(i) = subplot(2, 3, i);
end


% Plot initial time-series data with maximum values highlighted
max_indices = zeros(6, 1);
for i = 1:6
    plot(ax(i), filtered_MVIC{i});
    hold(ax(i), 'on');
    title(ax(i),MVIC_col_name(i))
    [max_val, max_idx] = max(filtered_MVIC{i});
    plot(ax(i), max_idx, max_val, 'ro');
    max_indices(i) = max_idx;
    yline(ax(i), max_val,'b:');
    ylim([0,max_val*1.1])
end

% Add button for storing selected values
uicontrol('Style', 'pushbutton', 'String', 'Store Values', 'Units','normalized',...
    'Position', [0.02 0.02 0.2 0.05], 'Callback', @store_values);

uicontrol(gcf,'Style','pushbutton','String','Done Selecting','Units','normalized',...
        'Position',[0.5 0.02 0.2 0.05],'Callback',@code_resume);
    
    % Define callback function for interactive selection
function select_callback(src,~, i)
    [x, y] = ginput(1);
    [~, idx] = min(abs(filtered_MVIC{i} - y));
    delete(findobj(ax(i), 'Type', 'line', 'Color', 'r', 'Marker', '*'));
    delete(findobj(ax(i), 'Type', 'ConstantLine', 'Color', 'g'));
    plot(ax(i),  idx, filtered_MVIC{i}(idx), 'r*');
    yline(ax(i), filtered_MVIC{i}(idx), 'g');
    max_indices(i) = idx;
end

% Add interactive selection to subplots
for i = 1:6
    set(ax(i), 'ButtonDownFcn', {@select_callback, i});
end
%

    % Define callback function for storing selected values
function store_values(src, event)
    selected_values = zeros(1, 6);
    for i = 1:6
        selected_values(i) = filtered_MVIC{i}(max_indices(i));
    end
    
    disp(selected_values);
    assignin('base', 'max_MVIC', selected_values);
    disp('Selected values saved as max_MVIC in workspace!');
    savefig(fig,'EMG_filter_MVIC_validation.fig')
    disp('Figures saved as EMG_filter_MVIC_validatio.fig!');
end

    function code_resume(src,event)
        uiresume(gcbf);
        close all;
    end

end