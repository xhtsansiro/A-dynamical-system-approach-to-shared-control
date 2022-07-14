%% plot the effect region of the VSDS background streamlines,
% based on the variance 

load("data_1.mat"); 

sigmascale = 1;
x_cen = (x_rec(:,1:end-1)+x_rec(:,2:end))/2;   % the center of the springs
x_len = vecnorm(x_rec(:,1:end-1)-x_rec(:,2:end)); % the length

weights_cal = @(x) omega_d(x, x_cen, x_len, sigmascale);
limits =  [-0.2, 0.5, 0, 0.5];
figure(1);
hold on
plot_colormap(limits, weights_cal);
