function [velocity_new] = velocity_limit(velocity)
% in linear dynamics, the far away from the equilbrium, the higher the 
% velocity is, which could be unfeasible. limit the velocity in 0.15m/s

    if norm(velocity)  > 0.20  % the max allow speed is 0.20 m/s
        velocity_new = 0.20 * velocity / norm(velocity);
    else
        velocity_new = velocity;
    end

end

