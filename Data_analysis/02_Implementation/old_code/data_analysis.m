%% This script performs data analysis on the simulation trial.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% The haptic device motion, real vs desired (sampling points)
hd_des = load('scale_down.txt');  % scaled traj from robot reference traj
hd_real = load('hd_real.txt');  % the real traj of haptic device
figure();
plot(hd_des(:,1), hd_des(:,2), 'b', 'Linewidth', 2);
hold on;
plot(hd_real(:,2), hd_real(:,3), 'r', 'Linewidth', 2);
xlabel('$y$', 'Interpreter','latex','FontSize',20);
ylabel('$z$', 'Interpreter','latex','FontSize',20);
grid on;
title(['Trajectory of the haptic device'], 'Interpreter','latex', 'FontSize',20);
legend({'$desired$', '$real$'}, 'Interpreter','latex', 'FontSize',20);
% xlim([-0.09, 0.08]); ylim ([-0.08, 0.10]);

%% The robot motion real vs desired 
robot_des = load('robot_desired.txt');   % desired one, map from the real movement of hd
robot_real = load('robot_real.txt');     % real movement of robot
figure();
plot(robot_des(:,2), robot_des(:,3), 'b', 'Linewidth', 2);
hold on;
plot(robot_real(:,2), robot_real(:,3), 'r', 'Linewidth', 2);
xlabel('$y$', 'Interpreter','latex','FontSize',20);
ylabel('$z$', 'Interpreter','latex','FontSize',20);
grid on;
title(['Trajectory of the Robot'], 'Interpreter','latex', 'FontSize',20);
legend({'$desired$', '$real$'}, 'Interpreter','latex', 'FontSize',20);
% xlim([-0.1, 0.5]); ylim ([0.1, 0.5]);

%% Plot the escape force which tries to get rid of the robot.
% vsds;
force = load('force.txt');   % passiv control force
figure();
tt = (0:length(force)-1) * 0.002;
plot(tt, force(:,1), 'b',  tt, force(:,2), 'r');
xlabel('$time:s$', 'Interpreter','latex','FontSize',20);
ylabel('$Force:N$', 'Interpreter','latex','FontSize',20);
grid on;
title(['Force acting on haptic device'], 'Interpreter','latex', 'FontSize',20);
legend({'$F_{y}$', '$F_{z}$'}, 'Interpreter','latex', 'FontSize',20);

% filter the trajectory
figure();
fs=500 ; dt=1/fs ;
[B,A]=butter(1,3/(fs/2)); 
force_filtered=filtfilt(B,A,force);
plot(tt, force_filtered(:,1), 'b',  tt, force_filtered(:,2), 'r');
xlabel('$time:s$', 'Interpreter','latex','FontSize',20);
ylabel('$Filtered Force:N$', 'Interpreter','latex','FontSize',20);
grid on;
title(['Filtered Force acting on haptic device'], 'Interpreter','latex', 'FontSize',20);
legend({'$F_{y}$', '$F_{z}$'}, 'Interpreter','latex', 'FontSize',20);


