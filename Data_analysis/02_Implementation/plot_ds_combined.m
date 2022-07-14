function [h,xd] = plot_ds_combined(fig, org_ds, new_ds, target, limits, par, varargin)
% this function plots the streamlines of ds + vsds(locally along reference traj.)
% par is the parameters, par.x_cen, par.x_len, par.w_th
    quality='medium';
    x_cen = par{1}; x_len = par{2}; w_th = par{3}; sigma_scale = par{4};

    if nargin > 6
        quality = varargin{1};
    end

    if strcmpi(quality,'high')
        nx=400;
        ny=400;
    elseif strcmpi(quality,'medium')
        nx=200;
        ny=200;
    else
        nx=50;
        ny=50;
    end

    axlim = limits;
    ax_x=linspace(axlim(1),axlim(2),nx); %computing the mesh points along each axis
    ax_y=linspace(axlim(3),axlim(4),ny); %computing the mesh points along each axis
    [x_tmp, y_tmp] = meshgrid(ax_x,ax_y); %meshing the input domain
    x=[x_tmp(:), y_tmp(:)]';  % laydown the matrix, x(:) reshape the matrix, 
    % xd = feval(ds, x-repmat(target,1,size(x,2))); % get the x_dot,
    %[xd, ~] = feval(ds, x); 
    xd = feval(org_ds, x);

    % check where the streamlines become VSDS
    omega_max = [];
    for i = 1:size(x,2)
        omega_max(i) = pos_check(x(:,i), x_cen, x_len, sigma_scale);
    end
% change those affected by 
    for i = 1:size(x,2)
        if omega_max(i) >= w_th
            xd(:,i) = feval(new_ds, x(:,i));
        end
    end

% the matrix containing dynamics are in a 2*n form, n are the samples
h = streamslice(x_tmp,y_tmp,reshape(xd(1,:),ny,nx),reshape(xd(2,:),ny,nx),4,'method','cubic');
% set(h,'LineWidth', 0.75)
set(h,'color',[0.0667  0.0667 0.0667]);

end

