function K_des = get_stiffness(x, var, stiff_type)
    % based on the variance of prediction, form the stiffness
    % x: position; var: variance; stiff_type: 
    
    k =  1/ var;  % scaling factor of stiffness,
    k = 1;  % 
    % the less variance is,  the bigger stiffness should be.
    if strcmp(stiff_type, 'variable')
        K_des = [50*(sin(2*x(2,1))+1)/2, 0; 0, 100]; % this sine function is for C U data, not robot
        
    elseif strcmp(stiff_type, 'constant')
        K_des = [25,0;0,100];
%         K_des = [100,0;0,200];
%         K_des = [500,0;0,1000];
%         K_des = k * [5,0;0,25];
    end
end


