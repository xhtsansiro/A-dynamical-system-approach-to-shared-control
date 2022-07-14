%% This script calculates the modulated params, scaling factor
% and rotation angle. 
% original dynamic is  t' = -5t， Vy = -5(y-0.4), Vz = -5(z-0.1)
%% Calculate velocity based on the sampled trajectory.
dt = 0.4; 
position = load("sampled.mat").sampled_point; 
% take the first 39 points, because the last point is inside the circle,
% where the originial dynamics works.
velocity =  diff(position) / dt;

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
    v_origin = [-5 * (position(i,1)-0.4); -5 * (position(i,2)-0.1)];
    v = [velocity(i,1); velocity(i,2)];
    theta = acos(dot (v_origin, v) / (norm(v_origin) * norm(v))); % in radian
    k = norm(v)/norm(v_origin) -1;
    p = [position(i,1), position(i,2), theta, k];
    input(i,:) =  p;
end

%% The input is the input for Gaussian process. 
% To calculate the velocity of the current states, it is important to know
% how many data points from the demonstration are nearby. Hence, a circle
% region centered at the current position will be made. 
start_pos = [0.2, 1]; 
% start_pos = [10, 10]; 
dtt = 0.002;
% start_pos = [-0.045, 0.14]; 
target = [0.4, 0.1];
% check how many points from the demonstrations are picked to formulate the
% dynamics at this position.
current = start_pos;
l_distance = 0.02; % this param should be tunable, 
sigma_f = 0.5; l = 1;
k_XX = kernel(input(:,1:2)); % avoid repeating of calculation.
mtx = inv(k_XX + 0.001*eye(39));  % inverse of the matrix 
trace = [current];
t = 0;
while t< 20
%while norm(current - target) > 0.03  % not inside the circle near equilibrium
    % check how many points from the demonstrations are picked to formulate the
    % dynamics at this position.
%     gp = [];
%     for i = 1: length(input)
%         k = 0;
%         if norm (current - [input(i,1), input(i,2)]) < l_distance
%             k = k+1;
%             gp(k,:) = input(i,:);
%         end
%     end
    % if gp does not have any element, which means the dynamics should be
    % the original one. 
    gp = input;
    v = -5 * (current- target);  % x' = -10(x-0.4)
    
    num = size(gp,1);   % how many rows
    if num ~= 0
        % gp has the information of nearby points, then formulate covariance
        % based on kernel functions. (miu, covariance)
        miu = zeros(length(gp),1);
        state = [gp(:,1:2);current];
        % optimize here, not necessary to calculate repeatedly. 
        for i = 1:size(state,1)-1
            k_xX(1,i) = 2 * exp(-(norm(state(i,:)- state(end,:)))^2/(2*0.05)); %l =1， sigma_f = 1
        end
        k_Xx = k_xX'; k_xx = 2 * exp(0); % sigma_f
       % k = kernel(state);
       % k_XX = k(1:num, 1:num); k_xX = k(num+1, 1:num);
       % k_Xx = k(1:num, num+1); k_xx = k(num+1, num+1);
        
        % apply GPR to find rotation angle and scaling factor.
        % A = inv(K_XX + 0.001*eye(num)); % noise 0.001  b/A = b*inv(A)
        %alpha = k_xX * inv( k_XX + 0.001*eye(num));
        alpha = k_xX * mtx;
        alpha = truncate(alpha);
        m_theta = alpha * gp(:,3);  % mean of the prediction, theta
        m_k = alpha * gp(:,4); % mean of prediction scaling factor
        % not allow m_k to be a very large value. 
        variance = k_xx - alpha * k_Xx;   % variance of the prediction
        
        % calculate the modulation 
        Mr = [cos(m_theta) -sin(m_theta); sin(m_theta) cos(m_theta)];
        M = (1+m_k) * Mr;
        v = (M * v')';
    end
    t = t + dtt; 
    current = current + dtt * v;
    % save the trace
    trace = [trace; current];
end

%% Plot the params 
figure(2)
scatter(trace(:,1), trace(:,2),8, 'g');
hold on
plot(trace(:,1), trace(:,2));
hold on
circle(0.4,0.1,0.02,2,'c');
hold on
plot(0.4,0.1, 'x', 'color','r', 'LineWidth',10)
grid on

hold on 
scatter(input(:,1), input(:,2))
