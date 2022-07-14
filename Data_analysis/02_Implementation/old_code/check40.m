%% Use sampling method to reduce the size, the interval between samples are 
% almost equal.

dt = 0.4; 
position = load("sampled.mat").sampled_point; 
% take the first 39 points, because the last point is inside the circle,
% where the originial dynamics works.
velocity =  diff(position) / dt;
start_pos = [0, 0.4];
% start_pos = [-0.045 , 0.15];   % starting from somewhere near the demonstration starti
% start_pos = [0.3, 0.3];
% start_pos = [2,2];
target = [0.4, 0.1];
% num = size(input,1);  

%% Plot velocity
tt = (1:39) * dt;
figure(1)
plot(tt, velocity(:,1),'b', tt, velocity(:,2), 'g')
grid on
legend('Vy', 'Vz')
xlabel('time: s')
ylabel('velocity: m/s')
title('Velocity curve')

%% Calculate rotation angle theta and scaling factor k 
input = [];
for i = 1:39
    % calculate original velocity, in the form of column vector
    v_origin = [-0.4 * (position(i,1)-0.4); -0.4 * (position(i,2)-0.1)];
    v = [velocity(i,1); velocity(i,2)];
    theta = acos(dot (v_origin, v) / (norm(v_origin) * norm(v))); % in radian
    k = norm(v)/norm(v_origin) -1;
    p = [position(i,1), position(i,2), theta, k];
    input(i,:) =  p;
end


%% Implement GPR process
current = start_pos; 

dt = 0.002;
miu = zeros(length(input),1);  % zero-mean
t = 0;
Br = 0.03; 
l_distance = 0.02;
l_distance_h = 0.04; 
idx = 0;

sigma_f = 5; l = 0.01;  % hyperparams for GP
variance = [];
trace = [current];
k_XX = sigma_f * exp( -l^-1 * pdist2(input(:,1:2), input(:,1:2)).^ 2 /2);
% k_XX = kernel(points(:,1:2));
% mtx = inv(k_XX + 0.0000001*eye(200));  % inverse of the matrix
k_XX = k_XX + 0.001*eye(39);   % no inverse of the matrix
R = chol(k_XX);  % do chelosky decomposition

%tmp_theta = R\(R'\points(:,3)) ;  % this is (k_XX + noise)^-1 * theta
%tmp_k = R\(R'\points(:,4)) ;
v_all = [];
k_Xx = zeros(1,39); k_xx = sigma_f; idx = 0;
while t < 40
%    idx = idx + 1;
%while norm (current-target) > 1e-5
%while norm(current - target) > 0.03
    
    idx = idx + 1;
%   state = [points(:,1:2);current];  % position information, used to calculate covariance matrixx
    
    k_Xx = sigma_f * exp( -l^-1 * pdist2 (current, input(:,1:2)).^2 /2); 
    
    alpha = (R\(R'\k_Xx'))';  % (k_XX ^-1 * k_xX)' = k_Xx * k_XX^-1
    alpha_t = truncate(alpha);   % truncate weights
%     alpha_t = alpha;
    theta_pred = alpha_t * input(:,3);
    k_pred = alpha_t * input(:,4);
    
    if k_pred + 1 < 1e-4
        k_pred = 0;
    elseif k_pred < -2
        k_pred = -2;
    elseif k_pred > 2
        k_pred = 2;
    else
        k_pred = 1* k_pred;
    end
    
    % calculate the prediction
   
%     theta_pred = k_Xx * tmp_theta;
%     k_pred = k_Xx * tmp_k;
%     tmp_variance = R\(R'\k_Xx');
    variance = k_xx - k_Xx * alpha';
    vv(idx) = variance;
    % alpha = k_Xx * mtx;
%     alpha = truncate(alpha);   % truncate weights
%     theta_pred = alpha * points(:,3);  % prediction of theta
%     scal_pred = alpha * points(:,4);  % prediction of scaling factor k
%     variance = k_xx - alpha * k_Xx'; 
    
    % calculate the velocity prediction based on theta and k
    Mr = [cos(theta_pred) -sin(theta_pred); sin(theta_pred) cos(theta_pred)];
    M = (1+k_pred) * Mr;
    v_o = -0.4 * (current- target);  % the originial linear dynamics
    v_n = (M * v_o')';   % modulated velocity
    v_all = [v_all; v_n];
    current = current + v_n * dt;

    t = t + dt;

    
        % if gp does not have any element, which means the dynamics should be
        % the original one. 
        % gp = input;
        
                
    % calculate the modulation 
%     Mr = [cos(m_theta) -sin(m_theta); sin(m_theta) cos(m_theta)];
%     M = (1+m_k) * Mr;
%     v = (M * v')';
%     
    % if there is no prediction at this point,
%     t = t + dt; 
%     current = current + dt * v;
%     % save the trace
    trace = [trace; current];
    
end


%% Plot the trajectory.
figure(2)
scatter(trace(:,1), trace(:,2),8, 'b');
hold on
plot(trace(:,1), trace(:,2));
hold on
circle(0.4,0.1,0.02,2,'c');
hold on
plot(0.4,0.1, 'x', 'color', 'r', 'LineWidth',10)
grid on

hold on 
scatter(input(:,1), input(:,2))
