%% This script uses butterworth filter to filter trajectory
% plot the filtered trajectory and the velocity.

traj = load('01_LfD/01_Data/trial3/data_mes.txt');
pos =[traj(:,2), traj(:,3)];  % first point put into the data of sampled points

% create butterworth filter to filter raw data.
fs=500 ; dt=1/fs ;
[B,A]=butter(1,3/(fs/2)); 
pos_filtered=filtfilt(B,A,pos); 

% plot the scatter plot of trajectory after filtering
figure()
plot(pos_filtered(:,1), pos_filtered(:,2), 'b', 'Linewidth',2);
hold on
% circle(0.4,0.1,0.02,2,'c');
% hold on
% plot(0.4,0.1, 'x', 'color','r', 'LineWidth',10);
scatter(0.4, 0.1, 150, [0 0 0],'d','Linewidth',2);
grid on 

xlabel('$x_y [m]$','Interpreter','LaTex','FontSize',20);
ylabel('$x_z [m]$','Interpreter','LaTex','FontSize',20);
% title(['Robot Trajectory of Demonstration'], 'Interpreter','LaTex','fontsize',20,'fontweight','normal');
xlim([-0.05,0.45]); ylim([0.05,0.45]);
legend({'$Trajectory$', '$Target$'}, 'Interpreter','latex', 'FontSize',10);
% title('Trajectory after filtering');

% create the velocity and filter it 
velocity = diff(pos_filtered)/dt;
velocity(end+1,:)=velocity(end,:) ;
t = (0: length(velocity)-1) * dt;
v_filtered = filtfilt(B,A,velocity); 

figure(2)
plot(t, v_filtered(:,1),'b')
hold on
plot(t, v_filtered(:,2),'r');
grid on
legend({'$v_y$', '$v_z$'}, 'Interpreter','latex', 'FontSize',16);
xlabel('$t:s$','Interpreter','LaTex','FontSize',20);
ylabel('$v:m/s$','Interpreter','LaTex','FontSize',20);
% title(['Velocity Profile before Segmentation'], 'Interpreter','LaTex','fontsize',20,'fontweight','normal');
%% Segmentataion to remove the part where the robot does not move at the 
% beginning and in the end. 
% start_pos= pos_filtered(1,:); %idx_end = pos_filtered(end,:);
flag_start = false;
for i = 1: length(v_filtered)
    if (norm(v_filtered(i,:)) > 1e-3 && flag_start == false)
        idx_start = i;
        flag_start = true;
    end
    
    % stop 
    if (flag_start == true  && norm(v_filtered(i,:)) < 1e-3)
        idx_end = i;
        break;
    end
    
end
tt = dt * (0:idx_end - idx_start);
pos_seg = pos_filtered(idx_start:idx_end,:);
v_seg = v_filtered(idx_start:idx_end, :);

figure(3)
plot(pos_seg(:,1), pos_seg(:,2),'b', 'Linewidth',2);
hold on
scatter(0.4, 0.1, 150, [0 0 0],'d','Linewidth',2);

xlabel('$x_y [m]$','Interpreter','LaTex','FontSize',20);
ylabel('$x_z [m]$','Interpreter','LaTex','FontSize',20);
xlim([-0.05,0.45]); ylim([0.05,0.45]);
legend({'$Trajectory$', '$Target$'}, 'Interpreter','latex', 'FontSize',15);

% title('Segemented Trajectory')
grid on

figure(4)
plot(tt, v_seg(:,1),'b')
hold on
grid on
plot(tt, v_seg(:,2), 'r');
legend({'$v_y$', '$v_z$'}, 'Interpreter','latex', 'FontSize',16);
xlabel('$t [s]$','Interpreter','LaTex','FontSize',20);
ylabel('$v [m/s]$','Interpreter','LaTex','FontSize',20);
% title(['Velocity Profile after Segmentation'], 'Interpreter','LaTex','fontsize',20,'fontweight','normal');

%% Use spline to downsample the data.

t_new = linspace(tt(1), tt(end), 200);  % downsampling to 200
pos_train(:,1) = spline(tt, pos_seg(:,1), t_new);
pos_train(:,2) = spline(tt, pos_seg(:,2), t_new);
v_train(:,1) = spline(tt, v_seg(:,1), t_new);
v_train(:,2) = spline(tt, v_seg(:,2), t_new);

%% Use sampled point to generate rotation angle and scaling factor.
% Calculate rotation angle theta and scaling factor k 
input = []; v_org = [];
for i = 1:200   % generate 200 points along the trajectory.
    % calculate original velocity, in the form of column vector
    v_origin = [-0.4 * (pos_train(i,1)-0.4); -0.4 * (pos_train(i,2)-0.1)];
    % linear dynamics is : x' = -0.4x;
    v = [v_train(i,1); v_train(i,2)];
    theta = acos(dot (v_origin, v) / (norm(v_origin) * norm(v))); % in radian
    % check it is clockwise or anti-clockwise rotation (from v_org to v_real) 
    if (v_origin(1) * v(2) - v_origin(2) * v(1) > 0)
        sgn = 1;
    else
        sgn = -1;
    end
    k = norm(v)/norm(v_origin) -1;
    v_org(:,end+1) = v_origin' ;
    p = [pos_train(i,1), pos_train(i,2), sgn*theta, k];
    input(i,:) =  p;
end
%%
