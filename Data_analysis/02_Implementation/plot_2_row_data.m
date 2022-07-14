function plot_2_row_data(data, leg1, leg2, til, namestr, lab1, lab2, dt)
    % plot
    tt = dt *(0:size(data,2)-1);
    figure();
    hold on;
    plot(tt, data(1,:), 'r', 'LineWidth',2);
    plot(tt, data(2,:), 'b', 'LineWidth',2);
    legend({leg1, leg2}, 'Interpreter','latex', 'FontSize',20);
 %   title([til, namestr], 'Interpreter','latex','FontSize',20);
    xlabel(lab1, 'Interpreter','latex', 'FontSize',20);
    ylabel(lab2, 'Interpreter','latex', 'FontSize',20);
    set(gca, 'Fontsize', 16);
    box on;
    grid on;
end