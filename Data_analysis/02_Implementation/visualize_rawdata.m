%% Load raw data from txt file.
A_soll = load('01_LfD/01_Data/trial2/data_comm.txt');
A_ist = load('01_LfD/01_Data/trial2/data_mes.txt');

B_soll = load('01_LfD/01_Data/trial3/data_comm.txt');
B_ist = load('01_LfD/01_Data/trial3/data_mes.txt');

C_soll = load('01_LfD/01_Data/trial4/data_comm.txt');
C_ist = load('01_LfD/01_Data/trial4/data_mes.txt');

% v_y = diff(A_soll(:,2))/0.002;
% v_z = diff(A_soll(:,3))/0.002;
% v = sqrt()
% t = (1:11258) * 0.002;
% plot(t, v_y);
% scatter(v_y, v_z,16,'b');

%% show the real trajectory of  raw data
figure(1)
scatter(A_ist(:,2), A_ist(:,3),8,'b');
hold on 
scatter(B_ist(:,2), B_ist(:,3),8, 'g');
hold on
scatter(C_ist(:,2), C_ist(:,3),8, 'k');
hold on 
circle(0.4,0.1,0.02,2,'c');
hold on
plot(0.4,0.1, 'x', 'color','r', 'LineWidth',10)
grid on

xlabel('y direction')
ylabel('z direction')
title('Ist trajectory')
%% show the target trajectory of the raw data
figure(2)
scatter(A_soll(:,2), A_soll(:,3),8,'b');
hold on 
scatter(B_soll(:,2), B_soll(:,3),8, 'g');
hold on
scatter(C_soll(:,2), C_soll(:,3),8, 'k');
grid on
xlabel('y direction')
ylabel('z direction')
title('Soll trajectory')

%% draw a circle around target point， which is (0.4, 0.1)
%circle(0.4,0.1,0.05,2,'c')



