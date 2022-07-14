function x_dot = linear_ds(x, A, att)
% linear dynamics, given the position, output the velocity field.
    x_dot =  A * (x-att);
end

