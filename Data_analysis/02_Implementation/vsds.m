function x_dot = vsds(x, A, x_rec, sigmascale, x0, th_begin)
% This script forms the dynamicas of vsds. 
% x_dot: the output of the vsds,
% x: position, col vector; x_rec: sampled points for VSDS; th_begin:
% threshold to avoid jerk behavior at the beginning. x0: starting position;
% A: the desired stiffness matrix, 
    K = size(x_rec,2) - 1;
    [N, M] = size(x);
    x_dot = zeros(N, M);
    fl = zeros(2, K);  
    x_cen = (x_rec(:,1:end-1)+x_rec(:,2:end))/2;   % the center of the springs
    x_len = vecnorm(x_rec(:,1:end-1)-x_rec(:,2:end)); % the length
    D = [0.8,0;0,0.8];
    
    omega = @(x) omega_t(x, x_cen, x_len, sigmascale); % calculate the weights
    
    for i = 1:M
        weights  = omega(x(:,i));  % return back a row vector
        alpha = startactivation(x(:,i), x0, th_begin);
        
        for k = 1:K
            % calculate the spring term from 1 to k-1 
            fl(:,k) = weights(k)*A(:,2*k-1:2*k)*(x(:,i)-x_rec(:,k+1));  % spring is used, starting from first point.
     
        end
        % at the end, each col of fl represents the spring force.
        fl_sum = sum(fl,2);
        % adding all elements of a row
        x_dot(:,i) = alpha * fl_sum;   

    end
    % plot_gaussian_along(x_cen, sigma);  % plot the gaussian along the trajectory
end


