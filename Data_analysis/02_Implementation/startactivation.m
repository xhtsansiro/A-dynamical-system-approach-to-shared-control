function alpha = startactivation(x,x0,th)
% act = 1.01 - exp(-(x-x0)'*(x-x0)/(0.15^2));


    b = 0.1; % 0.01 used before
    if norm(x-x0) > th   % th is shown as d in the paper 
        alpha = 1;
    else
        temp = asin(1-b);
        tempp = (temp)*(norm(x-x0)/th);
        alpha = sin(tempp) + b;
    end
end

