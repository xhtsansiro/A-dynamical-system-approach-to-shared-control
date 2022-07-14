function [alpha_new] = check(alpha)
%truncate the weights
    alpha_new = alpha;
    alpha_bar  = 0.010;  % 0.01 in paper  0.18: better result
    p = 0.010;  % 0.01 in paper  0.18: better result
    

        if alpha < 0
            sign = -1;
        else
            sign = 1;
        end
        alpha = sign * alpha(i);
        
        if alpha(i) < alpha_bar
            alpha_new(i) = 0;
        elseif alpha_bar <= alpha(i) <= alpha_bar + p
            alpha_new(i) = sign * 1/2*(1+ sin(pi*(alpha(i)- alpha_bar)/p - pi/2)) * alpha(i);
        else
            alpha_new(i) = sign * alpha(i);
        end
end