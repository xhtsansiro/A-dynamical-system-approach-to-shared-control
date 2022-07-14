%% load the points
tic
points = load('0.4velocity.mat').input; % velocity 0.4
att = [0.4; 0.1]; 
points_new = []; idx = 0;
for i = 1: size(points,1)
    idx = idx + 1;
    if norm(points(i,1:2)'-att) > 0.15    % 
        points_new(idx,:) = points(i,:);
    end
end

%% Form the gp dynamics, and plot streamlines
% hold on
att = [0.4; 0.1];
h_act = figure();
ds.limits = [-0.2, 0.7, -0.2, 0.7];
% ds.limits = [-0.3, 1.2, -0.4, 1.1]; % from -1 to 1 in both two dimensions
% att = [0; 0];  % the global attractor is at (0,0).
A  = - 0.4 * eye(2);
% myds = @(x) linear_dynamics(x, A_hat, att);  % linear dynamics
myds =  @(x) gp_ds (x, points_new, att);

scatter(points_new(:,1), points_new(:,2), 'b');
hold on; 
scatter(att(1),att(2), 150, [0 0 0],'d','Linewidth',2); 
hold on;
plot_ds_model(h_act, myds, [0;0], limits,'medium');

xlabel('$x_y [m]$','Interpreter','LaTex','FontSize',20);
ylabel('$x_z [m]$','Interpreter','LaTex','FontSize',20);
box on;
legend({'$selected\, points$', '$target$'}, 'Interpreter','latex', 'FontSize',15);

% f_VSDS is a function handle: f(x) 

%% plot the streamlines of linear DS
A  = - 0.4 * eye(2);
att = [0.4; 0.1]; 
h_act = figure();
lin_ds = @(x) l_ds(x, A, att);
ds.limits = [-0.2, 0.7, -0.2, 0.7];

[hatt_rob] = scatter(att(1),att(2), 150, [0 0 0],'d','Linewidth',2); hold on;
% [0 0 0] is the color RGB setting. 'd' represents the type of the marker.
[hds_rob1] = plot_ds_model(h_act, lin_ds, att, ds.limits,'medium'); hold on;

xlabel('$x_y [m]$','Interpreter','LaTex','FontSize',20);
ylabel('$x_z [m]$','Interpreter','LaTex','FontSize',20);
% title(['Streamlines of GP-based LMDS'], 'Interpreter','latex', 'FontSize',20);
set(gca,'fontsize',16,'LineWidth',1);
legend({'$target$'}, 'Interpreter','latex', 'FontSize',15);

box on;
axis(ds.limits);
 