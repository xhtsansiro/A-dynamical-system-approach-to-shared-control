%% This script performs Gaussian Process Regression based on the sampled points
points = load('0.4velocity.mat').input;

% start_pos = [-0.1, 0.7]; 
% start_pos = [-0.045 , 0.15];   % starting from somewhere near the demonstration starting point
 start_pos = [0, 0.2];
target = [0.4, 0.1];
num = size(points,1);   % how many rows

%% Implement GPR process
current = start_pos; 
traj = current;
dt = 0.002;
miu = zeros(length(points),1);  % zero-mean
t = 0;
Br = 0.03; 
l_distance = 0.08; % 0.02
l_distance_h = 0.04; 
idx = 0;
variance = [];
trace = [current];
while t < 20
% while norm (current-target) > 1e-5
% while norm(current - target) > 0.03
%     state = [points(:,1:2);current];  % position information, used to calculate covariance matrixx
%     % k = kernel(state);
%     for i = 1: size(state, 1)
%         
%     k_XX = k(1:num, 1:num); k_xX = k(num+1, 1:num);        
%     k_Xx = k(1:num, num+1); k_xx = k(num+1, num+1);
%     
%     alpha = k_xX / (k_XX + 0.001*eye(num));   % b/A = b*inv(A)
%     
%     % truncate weights
%     alpha = truncate(alpha);
%     theta_pred = alpha * points(:,3);  % prediction of theta
%     scal_pred = alpha * points(:,4);  % prediction of scaling factor k
%     variance = k_xx - alpha * k_Xx; 
%     
%     % calculate the velocity prediction based on theta and k
%     Mr = [cos(theta_pred) -sin(theta_pred); sin(theta_pred) cos(theta_pred)];
%     M = (1+scal_pred) * Mr;
%     v_o = -5 * (current- target);  % the originial linear dynamics
%     v_n = (M * v_o')';   % modulated velocity
%     current = current + v_n * dt;
%     traj= [traj; current];
%     t = t + dt;

    v = -0.4 * (current- target);  % x' = -10(x-0.4)

    if norm(current - target) <= Br
        m_theta = 0; m_k = 0;
    else
        gp = [];
        for i = 1: length(points)
            k = 0;
            if norm (current - [points(i,1), points(i,2)]) < l_distance
                k = k+1;
                gp(k,:) = points(i,:);
            end
        end
        % if gp does not have any element, which means the dynamics should be
        % the original one. 
        % gp = input;
  %      v = -5 * (current- target);  % x' = -10(x-0.4)
        idx = idx + 1;
    
        num = size(gp,1);   % how many rows
        if num ~= 0
            % gp has the information of nearby points, then formulate covariance
            % based on kernel functions. (miu, covariance)
            miu = zeros(length(gp),1);
            state = [gp(:,1:2);current];
            % optimize here, not necessary to calculate repeatedly. 
      
            k = kernel(state);
            k_XX = k(1:num, 1:num); k_xX = k(num+1, 1:num); 
            k_Xx = k(1:num, num+1); k_xx = k(num+1, num+1);
            % k_xX row vector, k_Xx col vector.
            
            % do the chelosky decomposition to calculate inverse
            R = chol(k_XX + 0.000001 * eye(num));
            alpha = (R\(R'\k_Xx))';  % (k_XX ^-1 * k_xX)' = k_Xx * k_XX^-1
            
            
            
            % apply GPR to find rotation angle and scaling factor.
            % A = inv(K_XX + 0.001*eye(num)); % noise 0.001  b/A = b*inv(A)
%           alpha = k_xX * inv( k_XX + 0.001*eye(num));
%           alpha = k_xX * mtx;
%           alpha = truncate(alpha);
            m_theta = alpha * gp(:,3);  % mean of the prediction, theta
            m_theta;
            m_k = alpha * gp(:,4); % mean of prediction scaling factor
           
            % limiting the scaling factor
            if m_k > 2
                m_k = 2;
            elseif m_k < -2
                m_k = -2;
            else
                m_k = 1 * m_k;
            end
            
            % not allow m_k to be a very large value. 
            variance(1,idx) = k_xx - alpha * k_Xx;   % variance of the prediction
        else
            m_theta = 0; m_k = 0;
            variance(1,idx) = 999; 
        end
%         % calculate the modulation 
%         Mr = [cos(m_theta) -sin(m_theta); sin(m_theta) cos(m_theta)];
%         M = (1+m_k) * Mr;
%         v = (M * v')';   
    end
                
%     gp = [];
%     for i = 1: length(points)
%         k = 0;
%         if norm (current - [points(i,1), points(i,2)]) < l_distance
%             k = k+1;
%             gp(k,:) = points(i,:);
%         end
%     end
%     % if gp does not have any element, which means the dynamics should be
%     % the original one. 
%     % gp = input;
%     v = -5 * (current- target);  % x' = -10(x-0.4)
%     idx = idx + 1;
%     
%     num = size(gp,1);   % how many rows
%     if num ~= 0
%         % gp has the information of nearby points, then formulate covariance
%         % based on kernel functions. (miu, covariance)
%         miu = zeros(length(gp),1);
%         state = [gp(:,1:2);current];
%         % optimize here, not necessary to calculate repeatedly. 
%       
%         k = kernel(state);
%         k_XX = k(1:num, 1:num); k_xX = k(num+1, 1:num);
%         k_Xx = k(1:num, num+1); k_xx = k(num+1, num+1);
%         
%         % apply GPR to find rotation angle and scaling factor.
%         % A = inv(K_XX + 0.001*eye(num)); % noise 0.001  b/A = b*inv(A)
%         alpha = k_xX * inv( k_XX + 0.001*eye(num));
% %         alpha = k_xX * mtx;
%         alpha = truncate(alpha);
%         m_theta = alpha * gp(:,3);  % mean of the prediction, theta
%         m_theta;
%         m_k = alpha * gp(:,4); % mean of prediction scaling factor
%         % not allow m_k to be a very large value. 
%         variance(1,idx) = k_xx - alpha * k_Xx;   % variance of the prediction
%     else
%         m_theta = 0; m_k = 0;
%         variance(1,idx) = 999; 
%     end
    % calculate the modulation 
    Mr = [cos(m_theta) -sin(m_theta); sin(m_theta) cos(m_theta)];
    M = (1+m_k) * Mr;
    v = (M * v')';
    
    % if there is no prediction at this point,
    t = t + dt; 
    current = current + dt * v;
    % save the trace
    trace = [trace; current];
    
end


%% Plot the trajectory.
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
scatter(points(:,1), points(:,2))




