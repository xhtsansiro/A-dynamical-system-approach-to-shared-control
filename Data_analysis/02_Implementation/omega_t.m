function omega = omega_t(x, x_cen, x_len, sigmascale)
    K = size(x_cen,2);
    omega = zeros(1,K);  % a row vector
    sigma = x_len * sigmascale;  % sigma in paper 
    for i = 1:K
        omega(i) = exp(-(1/(2*sigma(i)*sigma(i)))*(x-x_cen(:,i))'*(x-x_cen(:,i)));
%         omiga(i) = exp(-(1/(2*delta(i)*delta(i)))*(x-x_cen(:,i))'*(x-x_cen(:,i)));
%         omiga(i) = exp(-(1/(2*delta(i)))*(x-x_cen(:,i))'*(x-x_cen(:,i)));
    end
    
%     omiga
%     omiga = exp(4*omiga)
    omega_sum = sum(omega);
    omega = omega/omega_sum;

end