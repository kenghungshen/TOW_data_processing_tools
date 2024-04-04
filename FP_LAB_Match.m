%addpath('C:\Users\ks49639\Documents\MATLAB\ezc3d_matlab')

%%
%clc
%clear

%%
% load standard trial
load('C:\Users\ks49639\Box\Neurorehabilitation and Biomechanics Lab\R21 TOW\Data Processing Tools\FP_LAB_match_test\FP_origin_standard.mat') % data from s02 origin_restart02

% load static calibration from current trial
Static = ezc3dRead(static_file_name);
%% Extract data from c3d
% not using the middle marker as the position is not consistent
marker_header = string(Static.parameters.POINT.LABELS.DATA);
LtTread1 = find(marker_header=="LeftTreadmill1");
RtTread1 = find(marker_header=="RightTreadmill1");
LtTread3 = find(marker_header=="LeftTreadmill3");
RtTread3 = find(marker_header=="RightTreadmill3");

Static_LtTread1 = mean(squeeze(Static.data.points(:,LtTread1,:))');
Static_LtTread3 = mean(squeeze(Static.data.points(:,LtTread3,:))');
Static_RtTread1 = mean(squeeze(Static.data.points(:,RtTread1,:))');
Static_RtTread3 = mean(squeeze(Static.data.points(:,RtTread3,:))');

%% Method 1: offset

Static_treadmill_marker_all = [Static_LtTread1; Static_LtTread3; Static_RtTread1; Static_RtTread3];

static_minus_origin = (Static_treadmill_marker_all - Origin_treadmill_marker_all([1,3,4,6],:));
FP_origin_offset = mean(static_minus_origin);

disp('Calculate using method1: Treadmill marker offset')
disp(['FP origin in static trial is: ',num2str(1e-3*FP_origin_offset)]) % ezc3d target unit = mm, V3D unit = m, *1e-3 to report m
disp(['unit: m'])

figure(2)
subplot(1,2,1)

    % marker
scatter3(Static_treadmill_marker_all(1,1),Static_treadmill_marker_all(1,2),Static_treadmill_marker_all(1,3),'filled','b'); hold on
scatter3(Static_treadmill_marker_all(2,1),Static_treadmill_marker_all(2,2),Static_treadmill_marker_all(2,3),'filled','b'); hold on
scatter3(Static_treadmill_marker_all(3,1),Static_treadmill_marker_all(3,2),Static_treadmill_marker_all(3,3),'filled','r'); hold on
scatter3(Static_treadmill_marker_all(4,1),Static_treadmill_marker_all(4,2),Static_treadmill_marker_all(4,3),'filled','r'); hold on

scatter3(Origin_treadmill_marker_all(1,1),Origin_treadmill_marker_all(1,2),Origin_treadmill_marker_all(1,3),'b'); hold on
scatter3(Origin_treadmill_marker_all(3,1),Origin_treadmill_marker_all(3,2),Origin_treadmill_marker_all(3,3),'b'); hold on
scatter3(Origin_treadmill_marker_all(4,1),Origin_treadmill_marker_all(4,2),Origin_treadmill_marker_all(4,3),'r'); hold on
scatter3(Origin_treadmill_marker_all(6,1),Origin_treadmill_marker_all(6,2),Origin_treadmill_marker_all(6,3),'r'); hold on

    % offset
plot3([Origin_treadmill_marker_all(1,1);Static_LtTread1(1)],[Origin_treadmill_marker_all(1,2);Static_LtTread1(2)],[Origin_treadmill_marker_all(1,3);Static_LtTread1(3)],'k'); hold on
plot3([Origin_treadmill_marker_all(3,1);Static_LtTread3(1)],[Origin_treadmill_marker_all(3,2);Static_LtTread3(2)],[Origin_treadmill_marker_all(3,3);Static_LtTread3(3)],'k'); hold on
plot3([Origin_treadmill_marker_all(4,1);Static_RtTread1(1)],[Origin_treadmill_marker_all(4,2);Static_RtTread1(2)],[Origin_treadmill_marker_all(4,3);Static_RtTread1(3)],'k'); hold on
plot3([Origin_treadmill_marker_all(6,1);Static_RtTread3(1)],[Origin_treadmill_marker_all(6,2);Static_RtTread3(2)],[Origin_treadmill_marker_all(6,3);Static_RtTread3(3)],'k'); hold on

	% LAB CS
plot_x = [0,0,0;100,0,0];
plot_y = [0,0,0;0,100,0];
plot_z = [0,0,0;0,0,100];
scatter3(0,0,0,'b'); hold on
plot3(plot_x(:,1),plot_x(:,2),plot_x(:,3),'r'); hold on
plot3(plot_y(:,1),plot_y(:,2),plot_y(:,3),'g'); hold on
plot3(plot_z(:,1),plot_z(:,2),plot_z(:,3),'b'); hold on
	% FP CS
plot_x = [FP_origin_offset;FP_origin_offset+[100,0,0]];
plot_y = [FP_origin_offset;FP_origin_offset+[0,100,0]];
plot_z = [FP_origin_offset;FP_origin_offset+[0,0,100]];
plot3(plot_x(:,1),plot_x(:,2),plot_x(:,3),'r'); hold on
plot3(plot_y(:,1),plot_y(:,2),plot_y(:,3),'g'); hold on
plot3(plot_z(:,1),plot_z(:,2),plot_z(:,3),'b'); hold on
scatter3(FP_origin_offset(1),FP_origin_offset(2),FP_origin_offset(3),'filled','b'); hold on
    %LAB FP offset
plot_x = [0,0,0;FP_origin_offset];
plot3(plot_x(:,1),plot_x(:,2),plot_x(:,3),'k.'); hold on

xlim([-1000 1000])
ylim([-1000 1000])
zlim([-1000 1000])
title('Method 1: offset')
%% Method 2: Treadmill CS transformation

    % compute treadmill CS
TreadCS_x = Static_LtTread3 - Static_LtTread1;
TreadCS_x = TreadCS_x/norm(TreadCS_x);
TreadCS_y_temp = Static_RtTread1 - Static_LtTread1;
TreadCS_y_temp = TreadCS_y_temp/norm(TreadCS_y_temp);
TreadCS_z = cross(TreadCS_x,TreadCS_y_temp);
TreadCS_y = cross(TreadCS_z,TreadCS_x);
TreadCS_origin = Static_LtTread1;

FP_origin = [TreadCS_x', TreadCS_y', TreadCS_z']*lab_origin_TreadCS';
FP_origin = TreadCS_origin + FP_origin';
disp('-------------------------------------------------------------------------')
disp('Calculate using method 2: Treadmill frame transformation')
disp(['FP origin in static trial is: ',num2str(FP_origin*1e-3)]) % ezc3d target unit = mm, V3D unit = m, *1e-3 to report m
disp(['unit: m'])

figure(2)
subplot(1,2,2)
        % tread CS
plot_x = [TreadCS_origin;TreadCS_origin+TreadCS_x*100];
plot_y = [TreadCS_origin;TreadCS_origin+TreadCS_y*100];
plot_z = [TreadCS_origin;TreadCS_origin+TreadCS_z*100];
plot3(plot_x(:,1),plot_x(:,2),plot_x(:,3),'r'); hold on
plot3(plot_y(:,1),plot_y(:,2),plot_y(:,3),'g'); hold on
plot3(plot_z(:,1),plot_z(:,2),plot_z(:,3),'b'); hold on
        % LAB CS
plot_x = [0,0,0;100,0,0];
plot_y = [0,0,0;0,100,0];
plot_z = [0,0,0;0,0,100];
plot3(plot_x(:,1),plot_x(:,2),plot_x(:,3),'r:'); hold on
plot3(plot_y(:,1),plot_y(:,2),plot_y(:,3),'g:'); hold on
plot3(plot_z(:,1),plot_z(:,2),plot_z(:,3),'b:'); hold on
        % FP CS
plot_x = [FP_origin;FP_origin+[100,0,0]];
plot_y = [FP_origin;FP_origin+[0,100,0]];
plot_z = [FP_origin;FP_origin+[0,0,100]];
plot3(plot_x(:,1),plot_x(:,2),plot_x(:,3),'r'); hold on
plot3(plot_y(:,1),plot_y(:,2),plot_y(:,3),'g'); hold on
plot3(plot_z(:,1),plot_z(:,2),plot_z(:,3),'b'); hold on

        % tread marker
scatter3(Static_treadmill_marker_all(1,1),Static_treadmill_marker_all(1,2),Static_treadmill_marker_all(1,3),'filled','k'); hold on
scatter3(Static_treadmill_marker_all(2,1),Static_treadmill_marker_all(2,2),Static_treadmill_marker_all(2,3),'filled','b'); hold on
scatter3(Static_treadmill_marker_all(3,1),Static_treadmill_marker_all(3,2),Static_treadmill_marker_all(3,3),'filled','r'); hold on
scatter3(Static_treadmill_marker_all(4,1),Static_treadmill_marker_all(4,2),Static_treadmill_marker_all(4,3),'filled','r'); hold on

xlim([-1000 1000])
ylim([-1000 1000])
zlim([-1000 1000])
title('Method 2: CS transformation')
