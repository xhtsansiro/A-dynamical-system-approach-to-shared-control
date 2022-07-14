function omega_d = omega_d(x, x_cen, x_len, sigmascale)
    omega_d = [];
    K = size(x_cen,2);
%     omega = zeros(1,K);  % a row vector
    sigma = x_len * sigmascale;  % sigma in paper 
    

    for j = 1:size(x,2)
        omega = zeros(1,K);  % a row vector
        for i = 1:K
            omega(i) = exp(-(1/(2*sigma(i)*sigma(i)))*(x(:,j)-x_cen(:,i))'*(x(:,j)-x_cen(:,i)));
%           omiga(i) = exp(-(1/(2*delta(i)*delta(i)))*(x-x_cen(:,i))'*(x-x_cen(:,i)));
%           omiga(i) = exp(-(1/(2*delta(i)))*(x-x_cen(:,i))'*(x-x_cen(:,i)));
        end
        omega_d(end+1) = max(omega);  % dominant omega
    end
    
%     omiga
%     omiga = exp(4*omiga)
    
%     r = exp(-10*norm(x-att));
%     for i = 1:K
%         omiga(i) = (1-r)*omiga(i);
%     end
end
