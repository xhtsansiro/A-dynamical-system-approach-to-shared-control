function [h,xd] = plot_ds(fig, ds, target, limits, varargin)
% this function plots the streamlines of a certain ds
quality='medium';

if nargin > 4
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
xd = feval(ds, x);
% the matrix containing dynamics are in a 2*n form, n are the samples
h = streamslice(x_tmp,y_tmp,reshape(xd(1,:),ny,nx),reshape(xd(2,:),ny,nx),4,'method','cubic');
% set(h,'LineWidth', 0.75)
set(h,'color',[0.0667  0.0667 0.0667]);

end
