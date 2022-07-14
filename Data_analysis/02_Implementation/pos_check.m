function omega_d = pos_check(x, x_cen, x_len, sigmascale)
    % omega_d = [];
    K = size(x_cen,2);
%    omega = zeros(1,K);  % a row vector
    sigma = x_len * sigmascale;  % sigma in paper 
    
   % for j = 1:size(x,2)
    omega = zeros(1,K);  % a row vector
    for i = 1:K
%             omega(i) = exp(-(1/(2*sigma(i)*sigma(i)))*(x(:,j)-x_cen(:,i))'*(x(:,j)-x_cen(:,i)));
        omega(i) = exp(-(1/(2*sigma(i)*sigma(i)))*(x-x_cen(:,i))'*(x-x_cen(:,i)));
%           omiga(i) = exp(-(1/(2*delta(i)))*(x-x_cen(:,i))'*(x-x_cen(:,i)));
    end
%         omega_d(end+1) = max(omega);  % dominant omega
    omega_d = max(omega);  % dominant omega
  %  end
    
end
