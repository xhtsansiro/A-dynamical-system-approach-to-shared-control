%% This script is used to sample points to generate training samples for LMDS
% what important is the trajectory, instead of the time dependence. 

%% The sampled trajectory is composed of 2 parts, first part is the demonstration one, 
% when it is inside the circle near the target, it follows the given linear
% global dynamics, t' = -10t (t = x-x_goal)
% one use the trajectory of the 3rd trial. 

traj = load('01_LfD/01_Data/trial3/data_mes.txt');
sampled_point =[traj(1,2), traj(1,3)];  % first point put into the data of sampled points
previous = [traj(1,2), traj(1,3)];  
target = [0.4, 0.1];   % goal position
dt = 0.05;  % discrete interval
t = 0;

i = 2;
while i < length(traj)
    
    while i < length(traj) && norm([traj(i,2), traj(i,3)] - previous) < 0.02
        i = i + 1;
    end
   % if 
    if norm([traj(i,2), traj(i,3)] - target) < 0.02
        previous = [traj(i,2), traj(i,3)];
        sampled_point = [sampled_point; previous];
        break
%         cur = [traj(i,2), traj(i,3)];
%         sampled_point = [sampled_point; cur];
%         while t < 1
%             Vx = -10 * (cur(1) - 0.4);
%             Vy = -10 * (cur(2) - 0.1);
%             cur = cur + dt* [Vx, Vy];  % update
%             sampled_point = [sampled_point; cur];
%             t = t + dt; 
%            if cur ==  target
%                  break
%            end
%         end
%         break
    else
        previous = [traj(i,2), traj(i,3)];
        sampled_point = [sampled_point; previous];
    end

       
end

%% Plot
figure(1)
scatter(sampled_point(:,1), sampled_point(:,2),8,'b');
hold on 
plot(sampled_point(:,1), sampled_point(:,2));
hold on
circle(0.4,0.1,0.02,2,'c');
hold on
plot(0.4,0.1, 'x', 'color','r', 'LineWidth',10)
grid on

xlabel('y direction')
ylabel('z direction')
title('Sampled trajectory')

%% Save the sampled trajectory.
scatter(traj(:,2), traj(:,3),8,'b');
