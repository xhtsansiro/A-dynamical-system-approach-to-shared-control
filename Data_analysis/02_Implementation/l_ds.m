function x_dot = l_ds(x, A, att)
%
%   Detailed explanation goes here
     for i = 1: size(x,2)
        pos =  x(:,i);  % a column vector of current position 
        v_o = linear_ds(pos, A, att);  % the originial linear dynamics
        v_o = velocity_limit(v_o);  % limit the velocity. 
        x_dot(:,i) =v_o; % no limit on velocity.
     end
end

