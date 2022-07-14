%% Use spline to downsample to 200 points, use 200 points for regression.
% points = load('100points.mat').input; % 
points = load('0.4velocity.mat').input;
start_pos = [-0.0468, 0.1508];  % origin of the trajectory.
target = [0.4, 0.1];
Br = 0.15; 

% preprocessing, delete points close to the global stable point
points_new = []; cc = 0;
for i = 1: size(points,1)

    if norm(points(i,1:2)-target) > Br   % Inside Br, no points needed 
        points_new(end+1,:) = points(i,:);
%     else  % for sparsity sampling
%         cc = cc + 1;
%         if mod(cc, 10) == 0
%             points_new(end+1,:) = points(i,:);
%         end
    end
end
% points_new = points;

num = size(points_new,1);   % how many rows

%% Implement GPR process
current = start_pos; 
traj = current;
dt = 0.002;
% dt = 0.01;
miu = zeros(length(points_new),1);  % zero-mean
t = 0;
Br = 0.03; 
l_distance = 0.02;
l_distance_h = 0.04; 
idx = 0;

v_all =[]; vv = [];
% sigma_f = 2; l = 0.0015;  % hyperparams for GP  % (2, 0.004)
% sigma_f = 1; l = 0.0005;   % the new params used.
sigma_f = 1; l = 0.001;
variance = [];
trace = current;
k_XX = sigma_f * exp( -l^-1 * pdist2(points_new(:,1:2), points_new(:,1:2)).^ 2 /2);
% k_XX = kernel(points(:,1:2));
% mtx = inv(k_XX + 0.0000001*eye(200));  % inverse of the matrix
k_XX = k_XX + 0.01*eye(num);   % no inverse of the matrix
R = chol(k_XX);  % do chelosky decomposition

%tmp_theta = R\(R'\points(:,3)) ;  % this is (k_XX + noise)^-1 * theta
%tmp_k = R\(R'\points(:,4)) ;

k_Xx = zeros(1,num); k_xx = sigma_f; idx = 0;
% while t < 40
%    idx = idx + 1;
T = [];
tic
while norm (current-target) > 1e-2
%while norm(current - target) > 0.03
     
    idx = idx + 1;
%   state = [points(:,1:2);current];  % position information, used to calculate covariance matrixx
    tic   % start of the cycle
    k_Xx = sigma_f * exp( -l^-1 * pdist2 (current, points_new(:,1:2)).^2 /2); 
    
    alpha = (R\(R'\k_Xx'))';  % (k_XX ^-1 * k_xX)' = k_Xx * k_XX^-1
    alpha_t = truncate(alpha);   % truncate weights
%     alpha_t = alpha;
    theta_pred = alpha_t * points_new(:,3);
    k_pred = alpha_t * points_new(:,4);
        
    if abs(k_pred + 1) < 1e-2  % before was 1e-2 
        k_pred = 0;
    elseif k_pred < -1    % the scaling factor is not allowed to be smaller than -1.
        k_pred = 0;
    elseif k_pred > 2
        k_pred = 2;
    else
        k_pred = 1* k_pred;
    end

    % calculate the prediction
    variance = k_xx - k_Xx * alpha';
    vv(end+1) = variance;
    % alpha = k_Xx * mtx;
%     alpha = truncate(alpha);   % truncate weights
%     theta_pred = alpha * points(:,3);  % prediction of theta
%     scal_pred = alpha * points(:,4);  % prediction of scaling factor k
%     variance = k_xx - alpha * k_Xx'; 
    
    % calculate the velocity prediction based on theta and k
    Mr = [cos(theta_pred) -sin(theta_pred); sin(theta_pred) cos(theta_pred)];
    M = (1+k_pred) * Mr;
    v_o = -0.4 * (current- target);  % the originial linear dynamics
    v_o = velocity_limit(v_o);  % limit the velocity. 
    v_n = (M * v_o')';   % modulated velocity
    v_n = velocity_limit(v_n);  % limit the velocity. 
    
    % v_all = [v_all; v_n];
    v_all(:, end+1) = v_n;
    current = current + v_n * dt;
    T(end+1) = toc;
    % toc   % end of the cycle, 
%     if idx == 1220
%         ccc = 1;
%     end
    t = t + dt;
    trace(end+1,:) = current;
    % trace = [trace; current];
    
end
% toc
%% Plot the trajectory.
hold on
figure()
%scatter(trace(:,1), trace(:,2),8, 'b');
% hold on
plot(trace(:,1), trace(:,2),'LineWidth',2, 'Color', 'm');
% hold on
% circle(0.4,0.1,0.02,2,'c');
% hold on
% plot(0.4,0.1, 'x', 'color', 'r', 'LineWidth',10);
grid on
%%
hold on 
% scatter(points(:,1), points(:,2), 'g');
plot(points(:,1), points(:,2), 'LineWidth',2, 'Color', 'b')
scatter(target(1),target(2), 150, [0 0 0],'d','Linewidth',2); 
ylim([0.05, 0.45]); % the limits of y axis

hold on 
scatter(points_new(:,1), points_new(:,2), 'm');
legend('predicted trajectory', 'original demonstration', 'target', 'selected points');
xlabel('$y: m$','Interpreter','LaTex','FontSize',20);
ylabel('$z: m$','Interpreter','LaTex','FontSize',20);
title(['GP Predictive Trajectory'], 'Interpreter','latex', 'FontSize',20);


%% Plot velocity
load ('velocity.mat');
% create butterworth filter to filter raw data.
fs=500 ;
[B,A]=butter(1,3/(fs/2)); 
v_all=filtfilt(B,A,v_all'); 
tt = dt * (1:length(vv)); 

t_new = linspace(tt(1), tt(end), 200);  % downsampling to 200
v_t(:,1) = spline(tt, v_all(:,1), t_new);
v_t(:,2) = spline(tt, v_all(:,2), t_new);

figure()
plot (v_train(:,1));
hold on
plot (v_t(:,1));

figure()
plot (v_train(:,2));
hold on 
plot (v_t(:,2));
