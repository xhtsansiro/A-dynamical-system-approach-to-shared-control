%% This script is used to plot a evolve of the trajectory.
% generate a video of the trajectory.
%% synchronize the plot of the haptic device and the robot trajectory.
% original traj.
d_m = load("org_traj/vsds.txt");
h_r = load("org_traj/hd_real.txt");
r_r = load("org_traj/robot_real.txt");
limits_r = [-0.2, 0.7, -0.2, 0.7]; 
plot_title_r = 'Trajectory of Robot'; vid_title_r = 'video/robot_traj_1';
% robot 
cas = 1;
realtimePlot(r_r(:,2:3),d_m(:,1:2),limits_r, plot_title_r, vid_title_r, 1);
% haptic_device
limits_h = [-0.1, 0.1, -0.1, 0.1]; 
plot_title_h = 'Trajectory of Haptic Device'; vid_title_h = 'video/hd_traj_1';
realtimePlot(h_r(:,2:3),d_m(:,3:4),limits_h, plot_title_h, vid_title_h, 0);

%% start from pos (0.6755,0.2442), failed one
d_m = load("fail_traj_1/vsds.txt");
h_r = load("fail_traj_1/hd_real.txt");
r_r = load("fail_traj_1/robot_real.txt");
limits_r = [-0.2, 0.7, -0.2, 0.7]; 
plot_title_r = 'Trajectory of Robot'; vid_title_r = 'video/robot_traj_2';
% robot 
realtimePlot(r_r(1:4680,2:3),d_m(:,1:2),limits_r, plot_title_r, vid_title_r,2);
% haptic_device
limits_h = [-0.1, 0.1, -0.1, 0.1]; 
plot_title_h = 'Trajectory of Haptic Device'; vid_title_h = 'video/hd_traj_2';
% realtimePlot(h_r(1:4680,2:3),d_m(:,3:4),limits_h, plot_title_h, vid_title_h,0);

%% start from pos (0.6755,0.2442), successful one

d_m = load("suc_traj_1/vsds.txt");
h_r = load("suc_traj_1/hd_real.txt");
r_r = load("suc_traj_1/robot_real.txt");
limits_r = [-0.2, 0.7, -0.2, 0.7]; 
plot_title_r = 'Trajectory of Robot'; vid_title_r = 'video/robot_traj_3';
% robot 
realtimePlot(r_r(:,2:3),d_m(:,1:2),limits_r, plot_title_r, vid_title_r, 3);
% haptic_device
limits_h = [-0.05, 0.15, -0.05, 0.15]; 
plot_title_h = 'Trajectory of Haptic Device'; vid_title_h = 'video/hd_traj_3';
% realtimePlot(h_r(:,2:3),d_m(:,3:4),limits_h, plot_title_h, vid_title_h);

%% start again from (-0.0468, 0.1508), with a new obstacle, failed one
d_m = load("fail_traj_2/vsds.txt");
h_r = load("fail_traj_2/hd_real.txt");
r_r = load("fail_traj_2/robot_real.txt");
limits_r = [-0.2, 0.7, -0.2, 0.7]; 
plot_title_r = 'Trajectory of Robot'; vid_title_r = 'video/robot_traj_4';
% robot 
realtimePlot(r_r(1:6440,2:3),d_m(:,1:2),limits_r, plot_title_r, vid_title_r,4);
% haptic_device
limits_h = [-0.1, 0.1, -0.1, 0.1]; 
plot_title_h = 'Trajectory of Haptic Device'; vid_title_h = 'video/hd_traj_4';
% realtimePlot(h_r(1:6440,2:3),d_m(:,3:4),limits_h, plot_title_h, vid_title_h);

%% start again from (-0.0468, 0.1508), with a new obstacle, successful one
d_m = load("suc_traj_2/vsds.txt");
h_r = load("suc_traj_2/hd_real.txt");
r_r = load("suc_traj_2/robot_real.txt");
limits_r = [-0.2, 0.7, -0.2, 0.7]; 
plot_title_r = 'Trajectory of Robot'; vid_title_r = 'video/robot_traj_5';
% robot 
realtimePlot(r_r(:,2:3),d_m(:,1:2),limits_r, plot_title_r, vid_title_r, 5);
% haptic_device
limits_h = [-0.1, 0.1, -0.1, 0.1]; 
plot_title_h = 'Trajectory of Haptic Device'; vid_title_h = 'video/hd_traj_5';
% realtimePlot(h_r(:,2:3),d_m(:,3:4),limits_h, plot_title_h, vid_title_h);

%% visualize the escpae and demonstration
% first the escape
% start from pos (0.6755,0.2442), failed one
d_m = load("fail_traj_1/vsds.txt");
h_r = load("fail_traj_1/hd_real.txt");
r_r = load("escape/robot_real.txt");
limits_r = [-0.2, 0.7, -0.2, 0.7]; 
plot_title_r = 'Trajectory of Robot'; vid_title_r = 'video/robot_traj_6';
% robot 
realtimePlot(r_r(:,2:3),d_m(:,1:2),limits_r, plot_title_r, vid_title_r,2);

%% visualize the new demonstration
d_m = load("fail_traj_1/vsds.txt");
r_r = load("escape/filter_position.txt");
mov_back = load("escape/mov_back.txt");
limits_r = [-0.2, 0.7, -0.2, 0.7]; 
plot_title_r = 'New Demonstrations'; vid_title_r = 'video/robot_traj_7';
% robot 
realtimePlot(r_r(:,1:2),mov_back(:,2:3),limits_r, plot_title_r, vid_title_r,0);

%% visualize the move back
d_m = load("fail_traj_1/vsds.txt");
mov_back = load("escape/mov_back.txt");
limits_r = [-0.2, 0.7, -0.2, 0.7]; 
plot_title_r = 'Trajectory of Robot'; vid_title_r = 'video/robot_traj_8';
% robot 
realtimePlot(mov_back(:,2:3),d_m(:,1:2),limits_r, plot_title_r, vid_title_r,0);
