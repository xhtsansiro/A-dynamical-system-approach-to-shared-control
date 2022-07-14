%% this script defines the function to calculate kernel, 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this function can be discarded by using pdist2 to calculate the distance,
% each row represents one sample.

function K = kernel(x)
%KERNEL form the covariance matrix between X1 and X2
l = length(x);
sigma_f = 2.0; ll = 0.05;
for i =1:1:l
   for j = 1:1:l
       K(i,j) = sigma_f * exp(-(norm(x(i,:)- x(j,:)))^2/(2*ll)); %l =2， sigma_f = 0.05
       % l =0.1, sigma_f =1.0 from the github
   end
end

end