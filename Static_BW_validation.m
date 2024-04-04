% addpath('C:\Users\ks49639\Documents\MATLAB\ezc3d_matlab')
% 
% %%
% clc
% clear

%% Get sum Fz (BW) from static trial

% static_file_name = ["Static.c3d"];

% weight_scale_data = 69.1;
% weight_scale_session = ["Base"]

BW_validation(static_file_name,weight_scale_data)

function BW_validation(static_file_name,weight_scale_data)
%static_file_name = convertStringsToChars(static_file_name);
static  = ezc3dRead(convertStringsToChars(static_file_name)); 

    %find signal location handle in the matrix
    static_analog_header = string(static.parameters.ANALOG.LABELS.DATA);

    F1Z = find(static_analog_header=="Force.Fz1");
    F2Z = find(static_analog_header=="Force.Fz2");
    % get signal
    static_F1Z_sig = static.data.analogs(:,F1Z);
    static_F2Z_sig = static.data.analogs(:,F2Z);
    
    
    static_FZ_sig = static_F1Z_sig + static_F2Z_sig;
    
    static_FZ_sig = static_FZ_sig*-1;
    mean_FZ = mean(static_FZ_sig)/9.81;
    
    figure(3)
    plot(static_FZ_sig./9.81,'k--o'); hold on
    yline(mean_FZ,'r');
    yline(weight_scale_data,'b');
    legend('static trial Fz','mean static Fz','weight scale data')
    
    disp(' ')
    disp(['static trial Fz = ',num2str(mean_FZ), 'kg'])
    disp(['weight scale data = ', num2str(weight_scale_data), 'kg'])
    disp(['static trial Fz - weight scale data = ',num2str(mean_FZ - weight_scale_data), 'kg'])
    disp(' ')
end