%% Plot the reference traj. from three different position, 
% background is the streamline of GP-LMDS

% plot the trajectory of the robot, 
r1 = load("pos_1/robot_real.txt");
r2 = load("pos_2/robot_real.txt");
r3 = load("pos_3/robot_real.txt");
figure(1)
plot(r1(:,2), r1(:,3),'g', 'LineWidth', 2); 
hold on
plot(r2(1149:end,2), r2(1149:end,3),'b', 'LineWidth', 2);
hold on
plot(r3(200:end,2), r3(200:end,3),'r', 'LineWidth', 2);

points = load('0.4velocity.mat').input; % velocity 0.4
att = [0.4; 0.1]; 
points_new = []; idx = 0;
for i = 1: size(points,1)
    idx = idx + 1;
    if norm(points(i,1:2)'-att) > 0.15    % 
        points_new(idx,:) = points(i,:);
    end
end

att = [0.4; 0.1];
h_act = figure(1);
ds.limits = [-0.2, 0.7, -0.2, 0.7];
% ds.limits = [-0.3, 1.2, -0.4, 1.1]; % from -1 to 1 in both two dimensions
% att = [0; 0];  % the global attractor is at (0,0).
A  = - 0.4 * eye(2);
myds =  @(x) gp_ds (x, points_new, att);
[hatt_rob] = scatter(att(1),att(2), 150, [0 0 0],'d','Linewidth',2); hold on;
% [0 0 0] is the color RGB setting. 'd' represents the type of the marker.
[hds_rob1] = plot_ds(h_act, myds, att, ds.limits,'medium'); hold on;
xlabel('$x_y [m]$','Interpreter','LaTex','FontSize',20);
ylabel('$x_z [m]$','Interpreter','LaTex','FontSize',20);
% title(['Streamlines of GP-based LMDS'], 'Interpreter','latex', 'FontSize',20);
set(gca,'fontsize',16,'LineWidth',1);
box on;
axis(ds.limits);

legend({'$Trajectory_1$', '$Trajectory_2$',  '$Trajectory_3$',  '$Target$' }, 'Interpreter','Latex', 'FontSize',12);

%% Plot the stiffness along the three trajectories
%%% plot the eploside reflecting stiffness
% hold on, and do the plots with other trajectory, change the number to 2,1
pos = load('pos_1/pos_field.txt'); data.pos = pos(1:end-1,:)';
vel = load('pos_1/vel_field.txt'); v_vel = [];
% start from demonstration points
for i = 1:2:size(vel,1)
    v_vel(end+1,:) = vel(i,:);
end
data.vel = v_vel' ;
stiff = load('pos_1/stiffness.txt'); stiffness = [];
for i = 1:2:size(stiff,1)
    stiffness(:,:,end+1) = stiff(i:i+1,1:2);
end
data.stiffness = stiffness(:,:,2:end);

figure(2);
scatter(0.4, 0.1, 150, [0 0 0],'d','Linewidth',2);  % target is (0.4, 0.1)
hold on
% hold on
for i= 1:8   % this is for pos_1
% for i= 1:length(data.pos) -1  % this is for pos_2 and pos_3
    plotGMM2([data.pos(1,i); data.pos(2,i)], -data.stiffness(1:2,1:2,i) * 0.000001,[0.4940 0.1840 0.5560], .8);
%     plotGMM2([s_aug.x_att(1,i);s_aug.x_att(2,i)], s_aug.K_regress(1:2,1:2,i) * 0.002,[ 0.8500    0.3250    0.0980], .8);

    hold on
%    quiver(s_aug(k).x_att(1,i),s_aug(k).x_att(2,i),-s_aug(k).DataF(1,i)*0.0035,-s_aug(k).DataF(2,i)*0.0035,0,'r','LineWidth',3)
    quiver(data.pos(1,i), data.pos(2,i), data.vel(1,i)*0.0035, data.vel(2,i)*0.0035, 0, 'r', 'LineWidth', 3);
    hold on
%     plot(s_aug(k).x_att(1,i),s_aug(k).x_att(2,i),'*k','LineWidth',3)
    plot(data.pos(1,i), data.pos(2,i),'*k', 'LineWidth',4);
    hold on
end

hold on
plotGMM2([data.pos(1,1); data.pos(2,1)], -data.stiffness(1:2,1:2,1) * 0.000001,[0.4940 0.1840 0.5560], .8);
hold on
grid on
plot(pos(:,1), pos(:,2), 'g', 'Linewidth', 4);  % g: pos1, b: pos2, r: pos3;
set(gca,'fontsize',20,'LineWidth',1);
xlabel('$x_y [m]$','Interpreter','LaTex','FontSize',35);
ylabel('$x_z [m]$','Interpreter','LaTex','FontSize',35);
% xlim([0.3, 0.7]);
box on;
ylim ([0.05, 0.5]);
legend({'$target$', '$stiffness$' }, 'Interpreter','LaTex', 'FontSize',20);
% legend({'$Demonstration$', '$Target Position$', '$Stiffness$' }, 'Interpreter','latex', 'FontSize',10);
%% Plot the escape force which tries to get rid of the robot, add EscpaeForce directory in the path.
% vsds;
force_1 = load('escape_force_1.txt');   % starting from intial position the same as demonstration
force_2 = load('escape_force_2.txt');    % starting from another intial position

% plot force separately
figure();
tt_1 = (0:length(force_1)-1) * 0.002;
plot(tt_1, force_1(:,1), 'b',  tt_1, force_1(:,2), 'r');
xlabel('$time:s$', 'Interpreter','latex','FontSize',20);
ylabel('$Force:N$', 'Interpreter','latex','FontSize',20);
grid on;
box on
title(['Force acting on haptic device'], 'Interpreter','latex', 'FontSize',20);
legend({'$F_{y}$', '$F_{z}$'}, 'Interpreter','latex', 'FontSize',20);

figure();
tt_2 = (0:length(force_2)-1) * 0.002;
plot(tt_2, force_2(:,1), 'b',  tt_2, force_2(:,2), 'r');
xlabel('$time:s$', 'Interpreter','latex','FontSize',20);
ylabel('$Force:N$', 'Interpreter','latex','FontSize',20);
grid on;
box on
title(['Force acting on haptic device'], 'Interpreter','latex', 'FontSize',20);
legend({'$F_{y}$', '$F_{z}$'}, 'Interpreter','latex', 'FontSize',20);

%
% filter the trajectory
fs=500 ; dt=1/fs ;
[B,A]=butter(1,3/(fs/2)); 
force_filtered_1 = filtfilt(B,A,force_1);
force_filtered_2 = filtfilt(B,A,force_2);

comb_force_f_1 = vecnorm(force_filtered_1,2,2);
comb_force_f_2 = vecnorm(force_filtered_2,2,2);

% calculate the force combination from both axes

figure();
scatter(tt_1(1,3833),comb_force_f_1(3833,1), 70, [0 0 0],'x','Linewidth',4); hold on;
scatter(tt_2(1,2638),comb_force_f_2(2638,1), 70, [0 0 0],'x','Linewidth',4); hold on;
% scatter(attractor(1),attractor(2), 150, [0 0 0],'d','Linewidth',2); hold on;

plot(tt_1, comb_force_f_1, 'b', 'LineWidth', 4);   %,  tt, force_filtered(:,2), 'r');
hold on;
plot(tt_2, comb_force_f_2, 'r', 'LineWidth', 4); 
set(gca,'fontsize',20,'LineWidth',1);
x1 = xlabel('$Time [s]$');
% xlabel('$time [s]$', 'Interpreter','latex','FontSize',20);
y1 = ylabel('$||Force|| [N]$');
set([x1 y1],'interpreter','Latex','fontsize',30);
% ylabel('$Force:N$', 'Interpreter','latex','FontSize',20);
grid on;
box on
% title(['Escaping Force'], 'Interpreter','latex', 'FontSize',20);
scatter(tt_1(1,3833),comb_force_f_1(3833,1), 70, [0 0 0],'x','Linewidth',4); hold on;
scatter(tt_2(1,2638),comb_force_f_2(2638,1), 70, [0 0 0],'x','Linewidth',4); hold on;
legend({'$EscapingPoint_1$', '$EscapingPoint_2$'}, 'Interpreter','LaTex', 'FontSize',15);
% legend({'$F_{y}$', '$F_{z}$'}, 'Interpreter','latex', 'FontSize',20);


%% plot the effect region of the VSDS background streamlines, 
% Effective tude of VSDS and other area with GP-LMDS
% based on the variance, which determines the varying omega_threshold

% intialize GP-LMDS
points = load('sampled_points.mat').points_new;
limits = [-0.2, 0.7, -0.2, 0.7];
attractor = [0.4; 0.1];
N_points = 20; 
% N_points = 9;
ifplot = 0; % plot out the figure
stiff_type = 'constant'; 
x0 = [-0.0468; 0.1508];  % starting from original location
% x0 = [0.6; 0.4];
gp_lmds =  @(x) gp_ds (x, points, attractor);
% intialize VSDS based on GP-LMDS, starting from point x0.
[A_hat, B_, b_hat, x_rec, f_dot, th_begin] = get_vsds_parameters(x0, attractor, N_points, ifplot, limits, gp_lmds, stiff_type);
sigma_scale = 1;
my_vsds = @(x) vsds(x, A_hat, x_rec, sigma_scale, x0, th_begin);

%target = [0.4; 0.1];
h_act = figure();
hold on;
quality='medium';
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
ax_x=linspace(axlim(1),axlim(2),nx); % computing the mesh points along each axis
ax_y=linspace(axlim(3),axlim(4),ny); % computing the mesh points along each axis
[x_tmp, y_tmp]=meshgrid(ax_x,ax_y);  % meshing the input domain
x=[x_tmp(:), y_tmp(:)]';
x_ = x;
% x_ = x-repmat(attractor,1,size(x,2));
for i = 1:size(x_,2)
    xd(:,i) = feval(my_vsds, x_(:,i));
end

% calculate the weights
omega_max = zeros(1,nx*ny);
x_cen = (x_rec(:,1:end-1)+x_rec(:,2:end))/2;   % the center of the springs
x_len = vecnorm(x_rec(:,1:end-1)-x_rec(:,2:end)); % the length

for i = 1:size(x,2)
    omega_max(i) = pos_check(x_(:,i), x_cen, x_len, sigma_scale);
%     if omega_max(i) > 2
%         omega_max(i) = 0;
%     end
end

% for omega_threshold smaller than 0.4, have to use gp-lmds to recalculate
% omega_threshold = 0.8;  % big value is 0.7
omega_threshold = 0.1;
for i = 1:size(x_,2)
    if omega_max(i) < omega_threshold
        xd(:,i) = feval(gp_lmds, x_(:,i));
    end
end
z_tmp = reshape(omega_max,nx,ny);
hcolor = pcolor(x_tmp,y_tmp,reshape(omega_max,nx,ny));
set(hcolor,'linestyle','none');

% load whiteCopperColorMap;
% colormap(flipud(cm));
% colormap(cubehelix(128,2.3,-0.61,1.97,0.75,[0.50,0.96],[0.20,0.90]));
% colormap(cubehelix([],1.53,-1.32,2.74,0.81,[0.27,1],[0.38,0.86])); % temporary good parameter

% colormap(cubehelix([],1.53,-1.32,2.74,0.81,[0.5,1],[0.38,0.86])); % temporary good parameter

% colormap(cubehelix([],1.53,-1.32,2.74,0.81,[0.27,1],[0.38,0.45])); % w_th=0.8
colormap(cubehelix([],1.53,-1.32,2.74,0.81,[0.37,1],[0.35,0.78]));  % w_th = 0.1

% cubehelix_view(h_act)

colorbar;
caxis([min(min(z_tmp)), max(max(z_tmp))]);
h = streamslice(x_tmp,y_tmp,reshape(xd(1,:),ny,nx),reshape(xd(2,:),ny,nx),3,'method','cubic');
set(h,'LineWidth', 1)
set(h,'color',[0.0667  0.0667 0.0667]);

xl = xlabel('$x_y [m]$');
yl = ylabel('$x_z [m]$ ');
set(gca,'fontsize',20,'LineWidth',1);
set([xl yl],'interpreter','Latex','fontsize',35);
xlim([axlim(1),axlim(2)])
ylim([axlim(3),axlim(4)])
box on;
% set(gcf, 'Renderer', 'opengl')joint_external_force_callback
axis(limits);
scatter(attractor(1),attractor(2), 150, [0 0 0],'d','Linewidth',2); 
plot(x_rec(1,:), x_rec(2,:), 'm', 'Linewidth', 4);

%% plot the success trajectory in the intial cases, for postion 1.
h_act = figure();
real_traj = load('org_success_traj_1/robot_real.txt');
hold on; plot(real_traj(:,2), real_traj(:,3), 'b', 'LineWidth', 4); % take the first 4847 entries, ignore the rest chaotic 

ref_traj = load('org_success_traj_1/vsds.txt');
hold on; plot(ref_traj(:,1), ref_traj(:,2), 'r', 'LineWidth', 4, 'LineStyle','--');
att = [0.4; 0.1]; 
scatter(att(1),att(2), 150, [0 0 0],'d','Linewidth',2); hold on;
% border points
border_1 = [0.0160397163099708	0.124644755528446;
0.0302144760970289	0.158988615145704
0.0441609031805275	0.192028993506922
0.0561718791186619	0.219397343692300
0.0709293963915129	0.244707595298385
0.0886362323979195	0.267857215751910
0.110851932004782	0.289902042900067
0.132588288746976	0.308193593056834
0.154901499358700	0.322575197479801
0.178811803957508	0.333033807243606
0.202477311262400	0.339171606971543
0.226394510514650	0.340717870250050
0.244526206652862	0.339966870206809
0.254850631206484	0.338044449206074
0.269363607848400	0.324806064779534
0.288545748023157	0.288666311250792
0.303805158879610	0.257649792501695
0.314329940208349	0.237093590682349
0.320451229585832	0.206650178595676
0.323985642571862	0.167578027839024
0.328172577401719	0.129885283801361];
border_2 = [-0.109670116309971	0.176539244471554
-0.0950388760970289	0.211975384854296
-0.0795081031805275	0.248615006493078
-0.0600425991186619	0.290038656307700
-0.0355157963915129	0.329356404701615
-0.00620523239791949	0.365330784248090
0.0256984679952180	0.395943957099933
0.0629073112530243	0.424986406943166
0.105288500641300	0.449202802520199
0.151740196042492	0.466312192756395
0.201902688737600	0.475170393028457
0.253721489485350	0.473944129749951
0.309225793347138	0.459591129793191
0.365059368793516	0.421731550793926
0.394382392151600	0.378343935220466
0.404978251976843	0.351947688749208
0.428834841120390	0.311162207498305
0.447936059791651	0.262498409317651
0.455902770414168	0.218851821404324
0.459154357428138	0.182591972160976
0.463305422598281	0.145218716198639
];
% plot the wall
hold on
rectangle('Position', [0.295,0,0.01,0.235],'LineWidth', 2, 'LineStyle','-'); 
x = [0.295  0.305 0.305 0.295]; y=[0 0 0.235 0.235]; fill(x,y,[0.6350 0.0780 0.1840]);
hold on 
% scatter(border_1(:,1),border_1(:,2),'m', 'LineWidth',3.5);
plot(border_1(:,1),border_1(:,2),'m', 'LineWidth',4, 'LineStyle', ':');
hold on
rectangle('Position', [0.505,0,0.01,0.235],'LineWidth', 2, 'LineStyle','-'); 
x = [0.505  0.515 0.515 0.505]; y=[0 0 0.235 0.235]; fill(x,y,[0.6350 0.0780 0.1840]);
rectangle('Position', [0.295,0,0.21,0.01],'LineWidth', 2, 'LineStyle','-'); 
x = [0.295  0.505 0.505 0.295]; y=[0 0 0.01 0.01]; fill(x,y,[0.6350 0.0780 0.1840]);
% plot the border line of VSDS, generated based on the VSDS sampled
% attractors, calculate the velocity direction, then finds out its
% perpendicular direction vector k, then x +/- delta_l *k should be the
hold on
plot(border_2(:,1),border_2(:,2),'m', 'LineWidth',4, 'LineStyle', ':');
% plot ds streamlines in background 
% intialize GP-LMDS
points = load('sampled_points.mat').points_new; % GP dataset
limits = [-0.2, 0.7, -0.2, 0.7];
attractor = [0.4; 0.1];
N_points = 20; 
% N_points = 9;
ifplot = 0; % plot out the figure
stiff_type = 'constant'; 
x0 = [-0.0468; 0.1508];  % starting from original location
% x0 = [0.6; 0.4];
gp_lmds =  @(x) gp_ds (x, points, attractor);
% intialize VSDS based on GP-LMDS, starting from point x0.
[A_hat, B_, b_hat, x_rec, f_dot, th_begin] = get_vsds_parameters(x0, attractor, N_points, ifplot, limits, gp_lmds, stiff_type);
sigma_scale = 1;
my_vsds = @(x) vsds(x, A_hat, x_rec, sigma_scale, x0, th_begin);

x_cen = (x_rec(:,1:end-1)+x_rec(:,2:end))/2;   % the center of the springs
x_len = vecnorm(x_rec(:,1:end-1)-x_rec(:,2:end)); % the length
w_th = 0.3;
par = {x_cen, x_len, w_th, sigma_scale};
[hds_rob1,~] = plot_ds_combined(h_act, gp_lmds, my_vsds, att, limits, par, 'medium'); hold on;
set(gca,'fontsize',20,'LineWidth',1);
xlabel('$x_y [m]$','Interpreter','LaTex','FontSize',30);
ylabel('$x_z [m]$','Interpreter','LaTex','FontSize',30);
% title(['Streamlines of GP-based LMDS'], 'Interpreter','latex', 'FontSize',20);
% set(gca,'fontsize',16,'LineWidth',1);
box on;
axis(limits);

% legend({'$robot\,motion$', '$reference$', '$target$', '$wall$', '$border$'}, 'Interpreter','LaTex', 'FontSize',12);

%% plot the succeed trajectory in the intial cases, for postion 2.
h_act = figure();
real_traj = load('org_success_traj_2/robot_real.txt');
hold on; plot(real_traj(:,2), real_traj(:,3), 'b', 'LineWidth', 4); % take the first 4847 entries, ignore the rest chaotic 

ref_traj = load('org_success_traj_2/vsds.txt');
hold on; plot(ref_traj(:,1), ref_traj(:,2), 'r', 'LineWidth', 4, 'LineStyle','--');
att = [0.4; 0.1]; 
scatter(att(1),att(2), 150, [0 0 0],'d','Linewidth',2); hold on;
% border points
border_1 = [0.520392420562758	0.425683771148256
0.511285606967354	0.388249537411083
0.502225483027192	0.351004027916645
0.493246181537628	0.314094897303176
0.484461096651823	0.277211764873254
0.475825047360989	0.238913177708899
0.466394356777362	0.201452026943042
0.456913738946770	0.165064418045912
0.447898583685113	0.127812389558617];
border_2 = [0.442659579437242	0.444594228851744
0.433552393032646	0.407158462588917
0.424492516972809	0.369913972083355
0.415421818462372	0.332625102696824
0.406414903348177	0.294784235126746
0.398244952639011	0.258440822291101
0.388975643222638	0.221609973056958
0.379158261053230	0.183881581954088
0.370149416314887	0.146655610441383
];
% plot the wall
hold on
rectangle('Position', [0.295,0,0.01,0.235],'LineWidth', 2, 'LineStyle','-'); 
x = [0.295  0.305 0.305 0.295]; y=[0 0 0.235 0.235]; fill(x,y,[0.6350 0.0780 0.1840]);
hold on
plot(border_1(:,1),border_1(:,2),'m', 'LineWidth',4, 'LineStyle', ':');
hold on
rectangle('Position', [0.505,0,0.01,0.235],'LineWidth', 2, 'LineStyle','-'); 
x = [0.505  0.515 0.515 0.505]; y=[0 0 0.235 0.235]; fill(x,y,[0.6350 0.0780 0.1840]);
rectangle('Position', [0.295,0,0.21,0.01],'LineWidth', 2, 'LineStyle','-'); 
x = [0.295  0.505 0.505 0.295]; y=[0 0 0.01 0.01]; fill(x,y,[0.6350 0.0780 0.1840]);
hold on
plot(border_2(:,1),border_2(:,2),'m', 'LineWidth', 4, 'LineStyle', ':');

% plot ds streamlines in background 
% intialize GP-LMDS
points = load('sampled_points.mat').points_new; % GP dataset
limits = [-0.2, 0.7, -0.2, 0.7];
attractor = [0.4; 0.1];
%N_points = 20; 
N_points = 9;
ifplot = 0; % plot out the figure
stiff_type = 'constant'; 
x0 = [0.481526; 0.435139];  % starting from original location
% x0 = [0.6; 0.4];
gp_lmds =  @(x) gp_ds (x, points, attractor);
% intialize VSDS based on GP-LMDS, starting from point x0.
[A_hat, B_, b_hat, x_rec, f_dot, th_begin] = get_vsds_parameters(x0, attractor, N_points, ifplot, limits, gp_lmds, stiff_type);
sigma_scale = 1;
my_vsds = @(x) vsds(x, A_hat, x_rec, sigma_scale, x0, th_begin);

x_cen = (x_rec(:,1:end-1)+x_rec(:,2:end))/2;   % the center of the springs
x_len = vecnorm(x_rec(:,1:end-1)-x_rec(:,2:end)); % the length
w_th = 0.7;
par = {x_cen, x_len, w_th, sigma_scale};
[hds_rob1] = plot_ds_combined(h_act, gp_lmds, my_vsds, att, limits, par, 'medium'); hold on;
set(gca,'fontsize',20,'LineWidth',1);
xlabel('$x_y [m]$','Interpreter','LaTex','FontSize',30);
ylabel('$x_z [m]$','Interpreter','LaTex','FontSize',30);
% title(['Streamlines of GP-based LMDS'], 'Interpreter','latex', 'FontSize',20);

box on;
axis(limits);

legend({'$robot\,motion$', '$reference$', '$target$', '$wall$', '$border$'}, 'Interpreter','LaTex', 'FontSize',20, 'Location','northwest');


%% plot incremental learning. 
% plot the failed trajectory first
h_act = figure();
attractor = [0.4; 0.1];
fail_traj = load('failed_traj_case1/robot_real.txt');
hold on; plot(fail_traj(1:4680,2), fail_traj(1:4680,3), 'b', 'LineWidth', 4); % take the first 4847 entries, ignore the rest chaotic 

ref_traj = load('failed_traj_case1/vsds.txt');
hold on; plot(ref_traj(:,1), ref_traj(:,2), 'r', 'LineWidth', 4, 'LineStyle','--');

hold on; scatter(fail_traj(4680,2),fail_traj(4680,3), 150, [0 0 1],'x','Linewidth',2); 
scatter(attractor(1),attractor(2), 150, [0 0 0],'d','Linewidth',2); hold on;
border_1 = [0.685173782482615	0.170320491536146
0.651344527557225	0.156871992677011
0.617486930324934	0.143413550309795
0.583676907012381	0.129972541042012
0.549875620482338	0.116535029615008
0.516139113219779	0.103124623019266
0.482331848670525	0.0896845178485886
0.448555910131839	0.0762571447532369];
border_2 = [0.655620217517385	0.244661508463854
0.621793472442775	0.231214007322989
0.587933069675066	0.217754449690205
0.554123092987619	0.204313458957988
0.520324379517662	0.190876970384992
0.486584886780222	0.177465376980734
0.452778151329475	0.164025482151411
0.419004089868161	0.150598855246763];

% plot the wall
hold on
rectangle('Position', [0.295,0,0.01,0.235],'LineWidth', 2, 'LineStyle','-'); 
x = [0.295  0.305 0.305 0.295]; y=[0 0 0.235 0.235]; fill(x,y,[0.6350 0.0780 0.1840]);
hold on
plot(border_1(:,1),border_1(:,2),'m', 'LineWidth',4, 'LineStyle', ':');
hold on
rectangle('Position', [0.505,0,0.01,0.235],'LineWidth', 2, 'LineStyle','-'); 
x = [0.505  0.515 0.515 0.505]; y=[0 0 0.235 0.235]; fill(x,y,[0.6350 0.0780 0.1840]);
rectangle('Position', [0.295,0,0.21,0.01],'LineWidth', 2, 'LineStyle','-'); 
x = [0.295  0.505 0.505 0.295]; y=[0 0 0.01 0.01]; fill(x,y,[0.6350 0.0780 0.1840]);
plot(border_2(:,1),border_2(:,2),'m', 'LineWidth',4, 'LineStyle', ':');

points = load('sampled_points.mat').points_new; % GP dataset
limits = [-0.2, 0.7, -0.2, 0.7];
attractor = [0.4; 0.1];
%N_points = 20; 
N_points = 8;
ifplot = 0; % plot out the figure
stiff_type = 'constant'; 
x0 = [0.670397; 0.207491];  % starting from original location
% x0 = [0.6; 0.4];
gp_lmds =  @(x) gp_ds (x, points, attractor);
% intialize VSDS based on GP-LMDS, starting from point x0.
[A_hat, B_, b_hat, x_rec, f_dot, th_begin] = get_vsds_parameters(x0, attractor, N_points, ifplot, limits, gp_lmds, stiff_type);
sigma_scale = 1;
my_vsds = @(x) vsds(x, A_hat, x_rec, sigma_scale, x0, th_begin);

x_cen = (x_rec(:,1:end-1)+x_rec(:,2:end))/2;   % the center of the springs
x_len = vecnorm(x_rec(:,1:end-1)-x_rec(:,2:end)); % the length
w_th = 0.7;
par = {x_cen, x_len, w_th, sigma_scale};
[hds_rob1] = plot_ds_combined(h_act, gp_lmds, my_vsds, attractor, limits, par, 'medium'); hold on;
set(gca,'fontsize',20,'LineWidth',1);
xlabel('$x_y [m]$','Interpreter','LaTex','FontSize',30);
ylabel('$x_z [m]$','Interpreter','LaTex','FontSize',30);
% title(['Streamlines of GP-based LMDS'], 'Interpreter','latex', 'FontSize',20);

box on;
axis(limits);

hold on; scatter(fail_traj(4680,2),fail_traj(4680,3), 150, [0 0 1],'x','Linewidth',4); 

legend({'$robot\, motion$', '$reference$', '$stop\, point$', '$target$', '$wall$', '$border$'}, 'Interpreter','LaTex', 'FontSize',20);

%% escape from VSDS, starting from pos (0.6755,0.2442), Plot streamlines of relearnt DS
% after demonstration, plot the new streamlines
new_pos = load('new_points_case1/pos_train.txt');
new_vel = load('new_points_case1/vel_train.txt');
% new_pos = load('pos_train.txt');
% new_vel = load('vel_train.txt');
input = []; 
for i = 1:length(new_pos)   % generate 200 points along the trajectory.
    % calculate original velocity, in the form of column vector
    v_origin = [-0.4 * (new_pos(i,1)-0.4); -0.4 * (new_pos(i,2)-0.1)];
    % linear dynamics is : x' = -0.4x;
    v = [new_vel(i,1); new_vel(i,2)];
    theta = acos(dot (v_origin, v) / (norm(v_origin) * norm(v))); % in radian
    % check it is clockwise or anti-clockwise rotation (from v_org to v_real) 
    if (v_origin(1) * v(2) - v_origin(2) * v(1) > 0)
        sgn = 1;
    else
        sgn = -1;
    end
    k = norm(v)/norm(v_origin) -1;

    p = [new_pos(i,1), new_pos(i,2), sgn*theta, k];
    input(end+1,:) =  p;
end

att = [0.4; 0.1];
h_act = figure();
scatter(new_pos(:,1), new_pos(:,2), 'm'); hold on
ds.limits = [-0.2, 0.7, -0.2, 0.7];

A  = -0.4 * eye(2);
% myds = @(x) linear_dynamics(x, A_hat, att);  % linear dynamics
myds =  @(x) gp_ds (x, input, att);
[hatt_rob] = scatter(att(1),att(2), 150, [0 0 0],'d','Linewidth',2); hold on;
% [0 0 0] is the color RGB setting. 'd' represents the type of the marker.
[hds_rob1] = plot_ds(h_act, myds, att, ds.limits,'medium'); hold on;
xlabel('$x_y [m]$','Interpreter','LaTex','FontSize',20);
ylabel('$x_z [m]$','Interpreter','LaTex','FontSize',20);
hold on
legend({'selected points', 'target'}, 'Interpreter','none', 'FontSize',12);
% scatter(new_pos(:,1), new_pos(:,2), 'm');

%% plot the reference traj. and the real traj. of the robot. starting from pos (0.6755,0.2442);
% streamlines in the background, this case is the succeed one.
h_act = figure();
s_traj = load('success_traj_case1/robot_real.txt');
hold on; plot(s_traj(:,2), s_traj(:,3), 'b', 'LineWidth', 4); % take the first 4847 entries, ignore the rest chaotic 

ref_traj = load('success_traj_case1/vsds.txt');
hold on; plot(ref_traj(:,1), ref_traj(:,2), 'r', 'LineWidth', 4, 'LineStyle','--');
att = [0.4; 0.1];
hold on; 
scatter(att(1),att(2), 150, [0 0 0],'d','Linewidth',2); 
% border points
border_1 = [0.722384161887038	0.243324920975033
0.707157562202573	0.285904710570567
0.684539699564819	0.329829287053085
0.646514852288287	0.374369281832491
0.601768369736375	0.401404882198313
0.540932652961095	0.416947197819111
0.486426028486735	0.409538445141695
0.442451192679656	0.391004177551681
0.388057265150465	0.347971710103763
0.365571010467676	0.301547578579212
0.354695435848655	0.243687396556113
0.358708571269488	0.218084919638534
0.355914229194402	0.186080269819895
0.348696071648142	0.149127177261343
];

border_2 = [0.609015838112962	0.203985079024967
0.599274437797427	0.233357289429433
0.588648300435182	0.257684712946915
0.581029147711713	0.273812718167509
0.562409630263625	0.288043117801687
0.551493347038905	0.297412802180889
0.531227971513265	0.298215554858305
0.504520807320344	0.288303822448320
0.494378734849535	0.292332289896237
0.481824989532324	0.271798421420788
0.474000564151345	0.256582603443887
0.478157428730513	0.206597080361466
0.473687770805598	0.163071730180105
0.466329928351858	0.125414822738657
];

% plot the wall
hold on
rectangle('Position', [0.295,0,0.01,0.235],'LineWidth', 2, 'LineStyle','-'); 
x = [0.295  0.305 0.305 0.295]; y=[0 0 0.235 0.235]; fill(x,y,[0.6350 0.0780 0.1840]);
hold on
plot(border_1(:,1),border_1(:,2),'m', 'LineWidth',4, 'LineStyle', ':');
hold on
rectangle('Position', [0.505,0,0.01,0.235],'LineWidth', 2, 'LineStyle','-'); 
x = [0.505  0.515 0.515 0.505]; y=[0 0 0.235 0.235]; fill(x,y,[0.6350 0.0780 0.1840]);
rectangle('Position', [0.295,0,0.21,0.01],'LineWidth', 2, 'LineStyle','-'); 
x = [0.295  0.505 0.505 0.295]; y=[0 0 0.01 0.01]; fill(x,y,[0.6350 0.0780 0.1840]);
plot(border_2(:,1),border_2(:,2),'m', 'LineWidth',4, 'LineStyle', ':');
% plot streamlines
limits = [-0.2, 0.7, -0.2, 0.7];
N_points = 14;
ifplot = 0; % plot out the figure
stiff_type = 'constant'; 
x0 = [0.6657 ; 0.223655];  % starting from original location

gp_lmds = myds; % myds is generated in the previous section
% intialize VSDS based on GP-LMDS, starting from point x0.
[A_hat, B_, b_hat, x_rec, f_dot, th_begin] = get_vsds_parameters(x0, att, N_points, ifplot, limits, gp_lmds, stiff_type);
sigma_scale = 1;
my_vsds = @(x) vsds(x, A_hat, x_rec, sigma_scale, x0, th_begin);

x_cen = (x_rec(:,1:end-1)+x_rec(:,2:end))/2;   % the center of the springs
x_len = vecnorm(x_rec(:,1:end-1)-x_rec(:,2:end)); % the length
w_th = 0.3;
par = {x_cen, x_len, w_th, sigma_scale};
[hds_rob1] = plot_ds_combined(h_act, gp_lmds, my_vsds, att, limits, par, 'medium'); hold on;
set(gca,'fontsize',20,'LineWidth',1);
xlabel('$x_y [m]$','Interpreter','LaTex','FontSize',30);
ylabel('$x_z [m]$','Interpreter','LaTex','FontSize',30);
% title(['Streamlines of GP-based LMDS'], 'Interpreter','latex', 'FontSize',20);
box on;
axis(limits);

% legend({'$robot\,motion$', '$reference$', '$target$', '$wall$', '$border$'}, 'Interpreter','LaTex', 'FontSize',12, 'Location', 'northwest');

%% one new obstacle is added, traj failed before incremental learning.
% plot the failure case
h_act = figure();
fail_traj = load('failed_traj_case2/robot_real.txt');
hold on; plot(fail_traj(1:6440,2), fail_traj(1:6440,3), 'b', 'LineWidth', 4); % take the first 4847 entries, ignore the rest chaotic 

ref_traj = load('failed_traj_case2/vsds.txt');
hold on; plot(ref_traj(:,1), ref_traj(:,2), 'r', 'LineWidth', 4, 'LineStyle','--');

hold on; scatter(fail_traj(6440,2),fail_traj(6440,3), 150, [0 0 1],'x','Linewidth',4); 
scatter(att(1),att(2), 150, [0 0 0],'d','Linewidth',2); hold on;
% border points
% border points
border_1 = [0.0160397163099708	0.124644755528446;
0.0302144760970289	0.158988615145704
0.0441609031805275	0.192028993506922
0.0561718791186619	0.219397343692300
0.0709293963915129	0.244707595298385
0.0886362323979195	0.267857215751910
0.110851932004782	0.289902042900067
0.132588288746976	0.308193593056834
0.154901499358700	0.322575197479801
0.178811803957508	0.333033807243606
0.202477311262400	0.339171606971543
0.226394510514650	0.340717870250050
0.244526206652862	0.339966870206809
0.254850631206484	0.338044449206074
0.269363607848400	0.324806064779534
0.288545748023157	0.288666311250792
0.303805158879610	0.257649792501695
0.314329940208349	0.237093590682349
0.320451229585832	0.206650178595676
0.323985642571862	0.167578027839024
0.328172577401719	0.129885283801361];
border_2 = [-0.109670116309971	0.176539244471554
-0.0950388760970289	0.211975384854296
-0.0795081031805275	0.248615006493078
-0.0600425991186619	0.290038656307700
-0.0355157963915129	0.329356404701615
-0.00620523239791949	0.365330784248090
0.0256984679952180	0.395943957099933
0.0629073112530243	0.424986406943166
0.105288500641300	0.449202802520199
0.151740196042492	0.466312192756395
0.201902688737600	0.475170393028457
0.253721489485350	0.473944129749951
0.309225793347138	0.459591129793191
0.365059368793516	0.421731550793926
0.394382392151600	0.378343935220466
0.404978251976843	0.351947688749208
0.428834841120390	0.311162207498305
0.447936059791651	0.262498409317651
0.455902770414168	0.218851821404324
0.459154357428138	0.182591972160976
0.463305422598281	0.145218716198639
];

% plot the wall
hold on
rectangle('Position', [0.295,0,0.01,0.235],'LineWidth', 2, 'LineStyle','-'); 
x = [0.295  0.305 0.305 0.295]; y=[0 0 0.235 0.235]; fill(x,y,[0.6350 0.0780 0.1840]);
hold on
plot(border_1(:,1),border_1(:,2),'m', 'LineWidth',4, 'LineStyle', ':');
hold on
rectangle('Position', [0.505,0,0.01,0.235],'LineWidth', 2, 'LineStyle','-'); 
x = [0.505  0.515 0.515 0.505]; y=[0 0 0.235 0.235]; fill(x,y,[0.6350 0.0780 0.1840]);
rectangle('Position', [0.295,0,0.21,0.01],'LineWidth', 2, 'LineStyle','-'); 
x = [0.295  0.505 0.505 0.295]; y=[0 0 0.01 0.01]; fill(x,y,[0.6350 0.0780 0.1840]);
rectangle('Position', [0.19,0,0.02,0.4345],'LineWidth', 2, 'LineStyle','-'); 
x = [0.19  0.21 0.21 0.19]; y=[0 0 0.4345 0.4345]; fill(x,y,[0.6350 0.0780 0.1840]);
plot(border_2(:,1),border_2(:,2),'m', 'LineWidth',4, 'LineStyle', ':');

% plot ds streamlines in background 
limits = [-0.2, 0.7, -0.2, 0.7];
N_points = 20;
ifplot = 0; % plot out the figure
stiff_type = 'constant'; 

% starting from original location
x0 = [-0.0468; 0.1508];
gp_lmds = myds; % myds is generated in the previous section
% intialize VSDS based on GP-LMDS, starting from point x0.
[A_hat, B_, b_hat, x_rec, f_dot, th_begin] = get_vsds_parameters(x0, att, N_points, ifplot, limits, gp_lmds, stiff_type);
sigma_scale = 1;
my_vsds = @(x) vsds(x, A_hat, x_rec, sigma_scale, x0, th_begin);

x_cen = (x_rec(:,1:end-1)+x_rec(:,2:end))/2;   % the center of the springs
x_len = vecnorm(x_rec(:,1:end-1)-x_rec(:,2:end)); % the length
w_th = 0.3;
par = {x_cen, x_len, w_th, sigma_scale};
[hds_rob1] = plot_ds_combined(h_act, gp_lmds, my_vsds, att, limits, par, 'medium'); hold on;
set(gca,'fontsize',20,'LineWidth',1);
xlabel('$x_y [m]$','Interpreter','LaTex','FontSize',30);
ylabel('$x_z [m]$','Interpreter','LaTex','FontSize',30);
% title(['Streamlines of GP-based LMDS'], 'Interpreter','latex', 'FontSize',20);

box on;
axis(ds.limits);

hold on; scatter(fail_traj(6440,2),fail_traj(6440,3), 150, [0 0 1],'x','Linewidth',4); 
% legend({'$robot\, motion$', '$reference$', '$stop\, point$', '$target$', '$wall$', '$border$'}, 'Interpreter','LaTex', 'FontSize',12);

%% escape from VSDS, starting from the original pos (-0.0468; 0.1508);
% after demonstration, plot the new streamlines
new_pos = load('new_points_case2/pos_train.txt');
new_vel = load('new_points_case2/vel_train.txt');
input = []; 
for i = 1:length(new_pos)   % generate 200 points along the trajectory.
    % calculate original velocity, in the form of column vector
    v_origin = [-0.4 * (new_pos(i,1)-0.4); -0.4 * (new_pos(i,2)-0.1)];
    % linear dynamics is : x' = -0.4x;
    v = [new_vel(i,1); new_vel(i,2)];
    theta = acos(dot (v_origin, v) / (norm(v_origin) * norm(v))); % in radian
    % check it is clockwise or anti-clockwise rotation (from v_org to v_real) 
    if (v_origin(1) * v(2) - v_origin(2) * v(1) > 0)
        sgn = 1;
    else
        sgn = -1;
    end
    k = norm(v)/norm(v_origin) -1;

    p = [new_pos(i,1), new_pos(i,2), sgn*theta, k];
    input(end+1,:) =  p;
end

att = [0.4; 0.1];
h_act = figure(); 
scatter(new_pos(:,1), new_pos(:,2), 'b'); hold on;

ds.limits = [-0.2, 0.7, -0.2, 0.7];
% ds.limits = [-0.3, 1.2, -0.4, 1.1]; % from -1 to 1 in both two dimensions
% att = [0; 0];  % the global attractor is at (0,0).
A  = -0.4 * eye(2);
% myds = @(x) linear_dynamics(x, A_hat, att);  % linear dynamics
myds =  @(x) gp_ds (x, input, att);
[hatt_rob] = scatter(att(1),att(2), 150, [0 0 0],'d','Linewidth',2); hold on;
% [0 0 0] is the color RGB setting. 'd' represents the type of the marker.
[hds_rob1] = plot_ds(h_act, myds, att, ds.limits,'medium'); hold on;
xlabel('$x_y [m]$','Interpreter','LaTex','FontSize',20);
ylabel('$x_z [m]$','Interpreter','LaTex','FontSize',20);
hold on
box on
legend({'$selected\, points$', '$target$'}, 'Interpreter','Latex', 'FontSize',12);
%  scatter(new_pos(:,1), new_pos(:,2), 'm');

%% plot the traj which succeed to reach the target starting from original position. 
%  after updating the previous knowledge
h_act = figure();
s_traj = load('success_traj_case2/robot_real.txt');
hold on; plot(s_traj(:,2), s_traj(:,3), 'b', 'LineWidth', 4 ); % take the first 4847 entries, ignore the rest chaotic 

ref_traj = load('success_traj_case2/vsds.txt');
hold on; plot(ref_traj(:,1), ref_traj(:,2), 'r', 'LineWidth', 4, 'LineStyle','--');

hold on; 
scatter(att(1),att(2), 150, [0 0 0],'d','Linewidth',2); 

% hold on; scatter(s_traj(end,2),s_traj(end,3), 150, [0 0 1],'x','Linewidth',2); 
border_1 = [0.0214821857318995	0.108853403517477
0.0467771077116422	0.151582897186624
0.0676302570940142	0.233737351954124
0.0618305206491733	0.257329399232879
0.0588087151751709	0.274396536555574
0.0582629094617113	0.288227229523915
0.0768529663801525	0.315812839066942
0.0968628926376013	0.341495258018740
0.118337295510267	0.365034628109703
0.141344015418654	0.386213301174808
0.153456701045470	0.398844980131892
0.178560947099020	0.409188280359318
0.204372696146212	0.415747465379572
0.221315879961433	0.419724901914347
0.249727102700971	0.414013203882843
0.288291191042563	0.400212299574899
0.311827135232683	0.393844061232444
0.313484486964821	0.391068363938524
0.322273383336525	0.377593875117956
0.324219818444390	0.364885205620271
0.327434400215555	0.334122579099050
0.327145873579624	0.297800370888623
0.326429857670937	0.268285416335222
0.319228750255537	0.221497323437831
0.316883468540486	0.175108799654339
0.318301479220549	0.135323274981939
];
border_2 = [-0.114989585731899	0.192372596482523
-0.0997023077116422	0.215955102813376
-0.0896190570940142	0.234196648045876
-0.0981691206491733	0.256990600767121
-0.0953113151751709	0.317373463444426
-0.0738804894617113	0.378438770476085
-0.0485865663801525	0.415135160933058
-0.0202526926376013	0.450508741981260
0.0114227044897325	0.484069371890297
0.0465713845813459	0.515124698825192
0.0972272989545298	0.548639019868109
0.144855052900980	0.565597719640682
0.195415303853788	0.575496534620428
0.256194120038567	0.575877098085653
0.303722897299029	0.564626796117157
0.338652808957437	0.552079700425101
0.389250864767317	0.533863938767556
0.455263513035179	0.482221636061476
0.482122616663475	0.414538124882044
0.483556181555610	0.379442794379729
0.487427599784445	0.332647420900950
0.486996126420377	0.290879629111377
0.484340142329063	0.242510583664778
0.478959249744463	0.212214676562170
0.476782531459514	0.180791200345661
0.478134520779451	0.142630725018061
];

% plot the wall
hold on
rectangle('Position', [0.295,0,0.01,0.235],'LineWidth', 2, 'LineStyle','-'); 
x = [0.295  0.305 0.305 0.295]; y=[0 0 0.235 0.235]; fill(x,y,[0.6350 0.0780 0.1840]);
hold on
plot(border_1(:,1),border_1(:,2),'m', 'LineWidth',4, 'LineStyle', ':');
hold on

rectangle('Position', [0.505,0,0.01,0.235],'LineWidth', 2, 'LineStyle','-'); 
x = [0.505  0.515 0.515 0.505]; y=[0 0 0.235 0.235]; fill(x,y,[0.6350 0.0780 0.1840]);
rectangle('Position', [0.295,0,0.21,0.01],'LineWidth', 2, 'LineStyle','-'); 
x = [0.295  0.505 0.505 0.295]; y=[0 0 0.01 0.01]; fill(x,y,[0.6350 0.0780 0.1840]);
% a new wall
% rectangle('Position', [0.1875,0,0.025,0.4345],'LineWidth', 2, 'LineStyle','-'); 
% x = [0.1875  0.2125 0.2125 0.1875]; y=[0 0 0.4345 0.4345]; fill(x,y,'g');
rectangle('Position', [0.19,0,0.02,0.4345],'LineWidth', 2, 'LineStyle','-'); 
x = [0.19  0.21 0.21 0.19]; y=[0 0 0.4345 0.4345]; fill(x,y,[0.6350 0.0780 0.1840]);
plot(border_2(:,1),border_2(:,2),'m', 'LineWidth',4, 'LineStyle', ':');

% plot streamlines
x0 = [-0.0468; 0.1508];
N_points = 25; stiff_type = 'constant';
gp_lmds = myds; % myds is generated in the previous section
ifplot = 0;
% intialize VSDS based on GP-LMDS, starting from point x0.
[A_hat, B_, b_hat, x_rec, f_dot, th_begin] = get_vsds_parameters(x0, att, N_points, ifplot, ds.limits, gp_lmds, stiff_type);
sigma_scale = 1;
my_vsds = @(x) vsds(x, A_hat, x_rec, sigma_scale, x0, th_begin);

x_cen = (x_rec(:,1:end-1)+x_rec(:,2:end))/2;   % the center of the springs
x_len = vecnorm(x_rec(:,1:end-1)-x_rec(:,2:end)); % the length
w_th = 0.3;
par = {x_cen, x_len, w_th, sigma_scale};
[hds_rob1] = plot_ds_combined(h_act, gp_lmds, my_vsds, att, ds.limits, par, 'medium'); hold on;
set(gca,'fontsize',20,'LineWidth',1);
xlabel('$x_y [m]$','Interpreter','LaTex','FontSize',30);
ylabel('$x_z [m]$','Interpreter','LaTex','FontSize',30);
% title(['Streamlines of GP-based LMDS'], 'Interpreter','latex', 'FontSize',20);
box on;
axis(ds.limits);
% legend({'$robot\, motion$', '$reference$', '$target$', '$wall$','$border$'}, 'Interpreter','LaTex', 'FontSize',12);

%% plot an example of stiffness profile
var = 0:0.01:1; y = [];
k =1;
for i = 0:0.01:1
   y(k) = smooth_fall(i, 0.1, 0.8, 1100, 700);
   k = k+1;
end
figure();
plot(var, y, 'LineWidth', 4);
% set(gca,'fontsize',20,'LineWidth',10);
grid on
set(gca,'fontsize',25,'LineWidth',1);
xlabel('$\sigma^2$','Interpreter','LaTex','FontSize',35);
ylabel('$k [N/m]$','Interpreter','LaTex','FontSize',35);
% set(gca,'fontsize',30,'LineWidth',1);


%% plot the original demonstrations
points = load('0.4velocity.mat').input;
plot(points(:,1), points(:,2), 'LineWidth',2, 'Color', 'b');
hold on;
scatter(target(1),target(2), 150, [0 0 0],'d','Linewidth',2); 
ylim([0.05, 0.45]); % the limits of y axis
grid on
box on

xlabel('$x_y [m]$','Interpreter','LaTex','FontSize',20);
ylabel('$x_z [m]$','Interpreter','LaTex','FontSize',20);
% title(['GP Predictive Trajectory'], 'Interpreter','latex', 'FontSize',20);


%% plot the demonstrations used for incremental learning

points = load('new_points_case1/refine.txt');
plot(points(:,1), points(:,2), 'LineWidth',2, 'Color', 'b');
hold on;
target = [0.4;0.1];
scatter(target(1),target(2), 150, [0 0 0],'d','Linewidth',2); 
ylim([0.05, 0.45]); % the limits of y axis
grid on
box on

xlabel('$x_y [m]$','Interpreter','LaTex','FontSize',20);
ylabel('$x_z [m]$','Interpreter','LaTex','FontSize',20);

%% data prcoessing;
new_pos = load('new_points_case1/pos_train.txt');
new_vel = load('new_points_case1/vel_train.txt');

%simplify points
points_take = [];
for i = 120:length(new_pos)
    if mod(i,2) == 0
        points_take(end+1,:) = new_pos(i,:);
    end
end
concat = [new_pos(1:120,:); points_take];

input = []; 
for i = 1:length(new_pos)   % generate 200 points along the trajectory.
    % calculate original velocity, in the form of column vector
    v_origin = [-0.4 * (new_pos(i,1)-0.4); -0.4 * (new_pos(i,2)-0.1)];
    % linear dynamics is : x' = -0.4x;
    v = [new_vel(i,1); new_vel(i,2)];
    theta = acos(dot (v_origin, v) / (norm(v_origin) * norm(v))); % in radian
    % check it is clockwise or anti-clockwise rotation (from v_org to v_real) 
    if (v_origin(1) * v(2) - v_origin(2) * v(1) > 0)
        sgn = 1;
    else
        sgn = -1;
    end
    k = norm(v)/norm(v_origin) -1;

    p = [new_pos(i,1), new_pos(i,2), sgn*theta, k];
    input(end+1,:) =  p;
end

att = [0.4; 0.1];
h_act = figure();
scatter(concat(:,1), concat(:,2), 'b'); hold on
ds.limits = [-0.2, 0.7, -0.2, 0.7];

A  = -0.4 * eye(2);
% myds = @(x) linear_dynamics(x, A_hat, att);  % linear dynamics
myds =  @(x) gp_ds (x, input, att);
[hatt_rob] = scatter(att(1),att(2), 150, [0 0 0],'d','Linewidth',2); hold on;
% [0 0 0] is the color RGB setting. 'd' represents the type of the marker.
[hds_rob1] = plot_ds(h_act, myds, att, ds.limits,'medium'); hold on;
xlabel('$x_y [m]$','Interpreter','LaTex','FontSize',20);
ylabel('$x_z [m]$','Interpreter','LaTex','FontSize',20);
hold on
legend({'$selected\,points$', '$target$'}, 'Interpreter','latex', 'FontSize',12);
box on

















