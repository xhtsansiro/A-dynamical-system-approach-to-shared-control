clc

robot_real_1 = load('start from demonstration points/robot_real.txt');     % real movement of robot
figure(2);
hold on
plot(robot_real_1(:,2), robot_real_1(:,3), 'g', 'Linewidth', 4);
xlabel('$y$', 'Interpreter','latex','FontSize',20);
ylabel('$z$', 'Interpreter','latex','FontSize',20);
grid on;
title(['Trajectory of the Robot'], 'Interpreter','latex', 'FontSize',20);
% legend({'$Trajectory$'}, 'Interpreter','latex', 'FontSize',20);

%% from position 2
robot_real_2 = load('start from Pos1/robot_real.txt');     % real movement of robot
hold on

plot(robot_real_2(:,2), robot_real_2(:,3), 'r', 'Linewidth', 4);
% xlabel('$y$', 'Interpreter','latex','FontSize',20);
% ylabel('$z$', 'Interpreter','latex','FontSize',20);
% grid on;
% title(['Trajectory of the Robot'], 'Interpreter','latex', 'FontSize',20);
% legend({'$Trajectory 2$'}, 'Interpreter','latex', 'FontSize',20);

%% from position 3
robot_real_3 = load('start from Pos2/robot_real.txt');     % real movement of robot
% figure();
hold on
plot(robot_real_3(:,2), robot_real_3(:,3), 'b', 'Linewidth', 4);
% xlabel('$y$', 'Interpreter','latex','FontSize',20);
% ylabel('$z$', 'Interpreter','latex','FontSize',20);
% grid on;
% title(['Trajectory of the Robot'], 'Interpreter','latex', 'FontSize',20);
% legend({'$Target$','$Trajectory 1$', '$Trajectory 2$', '$Trajectory 3$'}, 'Interpreter','latex', 'FontSize',15);
legend({'$Trajectory 1$', '$Trajectory 2$', '$Trajectory 3$'}, 'Interpreter','latex', 'FontSize',15);