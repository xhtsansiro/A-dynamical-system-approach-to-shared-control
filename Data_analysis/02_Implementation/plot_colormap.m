function [h] = plot_colormap(limits,func)
%PLOT_COLORMAP Summary of this function goes here
    nx = 400; ny = 400;
    xlim  = limits;
    ax_x = linspace(xlim(1), xlim(2), nx);
    ax_y = linspace(xlim(3), xlim(4), ny);
    [x_tmp, y_tmp] = meshgrid(ax_x, ax_y);
    x = [x_tmp(:), y_tmp(:)]';
    xd = feval(func, x);
    z_tmp = reshape(xd, nx, ny);

    h = pcolor(x_tmp,y_tmp,z_tmp);
    set(h, 'linestyle', 'none');
    colormap(jet);
    colorbar;
end

