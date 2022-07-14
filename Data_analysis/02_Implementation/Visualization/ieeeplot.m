%% plot simulated stiffness
K_stiff = [];

for i = 1:10:size(x,2)
    K_stiff(:,end+1) = plot_stiffness_elipse(x(:,i), grad_ds, myds);
end
a = min(x(1,:));
b = max(x(1,:));

figure;
hold on;
plot(K_stiff(1,:));
plot(K_stiff(2,:));

figure
hold on;
xx = x(1,1:10:size(x,2));
lok = 3000;
plot(xx(1:lok),K_stiff(1,1:lok),'r','linewidth',2);
plot(xx(1:lok),K_stiff(2,1:lok),'b','linewidth',2);
% plot(xx(1:lok), 50*(sin(8*xx(1:lok)+1)+1)/2, 'r--','linewidth',2);
% plot(xx(1:lok), 100*(sin(6*xx(1:lok)+4)+2)/2, 'b--','linewidth',2);
% plot(xx(1:lok), 25*ones(1,lok), 'r--','linewidth',2);
% plot(xx(1:lok), 100*(sin(6*xx(1:lok)+4)+2)/2, 'b--','linewidth',2);
% plot(xx(1:lok), 50*(sin(8*xx(1:lok)+1)+1)/2, 'r--','linewidth',2);
plot(xx(1:lok), 100*ones(1,lok), 'b--','linewidth',2);
plot(xx(1:lok), 25*ones(1,lok), 'r--','linewidth',2);
set(gca,'fontsize',20,'LineWidth',1);
% title("Stiffness Profile", 'Interpreter','LaTex','FontSize',20)
xlabel('$x_y [m]$','Interpreter','LaTex','FontSize',24);
ylabel('$K [N/m]$','Interpreter','LaTex','FontSize',24);
box on;
legend('K_1','K_2','K_{des,1}','K_{des,2}')
% xlim([0 sizeofdata])
ylim([0 120])
xlim([a b])


%% plot SEDS and VSDS 
h_act = figure();
[hatt_rob] = scatter(att(1),att(2), 150, [0 0 0],'d','Linewidth',2); hold on;
[hds_rob1] = plot_ds_model(h_act, myds, [0;0], limits,'medium'); hold on;
xlabel('$x_1$','Interpreter','LaTex','FontSize',24);
ylabel('$x_2$','Interpreter','LaTex','FontSize',24);
set(gca,'fontsize',20);
box on;
axis(limits);
set(gca,'XTick',[], 'YTick', [])
plot(x_rec(1,:),x_rec(2,:),'r','linewidth',4);
plot(x_rec(1,1),x_rec(2,1),'b*','markersize',10,'Linewidth',2);

h_act = figure();
[hatt_rob] = scatter(att(1),att(2), 150, [0 0 0],'d','Linewidth',2); hold on;
[hds_rob1] = plot_ds_model(h_act, global_ds, [0;0], limits,'medium'); hold on;
xlabel('$x_1$','Interpreter','LaTex','FontSize',24);
ylabel('$x_2$','Interpreter','LaTex','FontSize',24);
set(gca,'fontsize',20);
box on;
axis(limits);
set(gca,'XTick',[], 'YTick', [])
plot(x_rec(1,:),x_rec(2,:),'r','linewidth',4);
plot(x_rec(1,1),x_rec(2,1),'b*','markersize',10,'Linewidth',2);

%% plot stiffness colormap in insertion task
t = -0.02:0.001:0.179;
m = (t+0.02)/(0.179+0.02);
k22 = zeros(1,length(t));
k11 = zeros(1,length(t));
t_max = 0.05;
t_min = 0.0;
k22_min = 100;
k22_diff = 600;
e_max = 0.05;
e_min = 0.0;
k11_min = 700;
k11_diff = 500;

h_act = figure();

box on;
axis(limits);
target = [0;0];
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
x_ = x-repmat(target,1,size(x,2));
omega_max = zeros(1,nx*ny);
delta = x_len * sigmascale;
for i = 1:size(x,2)
%     omega_max(i) = k11_min + k11_diff*smooth_transition_fall(x_(2,i) ,e_max,e_min);
    omega_max(i) =  k22_min + smooth_transition_rising(x_(2,i),t_max,t_min,k22_diff);
end
z_tmp = reshape(omega_max,nx,ny);
hcolor = pcolor(x_tmp,y_tmp,reshape(omega_max,nx,ny));
set(hcolor,'linestyle','none');
colormap(cubehelix([],1.53,-1.32,2.74,0.81,[0.27,1],[0.38,0.86])); % temporary good parameter
colorbar;
caxis([min(min(z_tmp)), max(max(z_tmp))]);
[hatt_rob] = scatter(att(1),att(2), 150, [0 0 0],'d','Linewidth',2); hold on;
[hds_rob1] = plot_ds_model(h_act, myds, [0;0], limits,'medium'); hold on;
set(gca,'fontsize',20);
xlabel('$x_y [m]$','Interpreter','LaTex','FontSize',24);
ylabel('$x_z [m]$','Interpreter','LaTex','FontSize',24);


%% plot crash force norm
tau_ext_mycrash = load('tau_ext_mycrash.txt');
tau_ext_sedscrash = load('tau_ext_sedscrash.txt');

fs=500 ; 
[B,A]=butter(1,1.2/(fs/2));
tau_ext_mycrash_filt=filtfilt(B,A,tau_ext_mycrash) ;
tau_ext_sedscrash_filt=filtfilt(B,A,tau_ext_sedscrash) ;

t = 0:0.002:11.998;

figure;
plot(t, vecnorm(tau_ext_mycrash_filt(1:6000,1:3)'),'linewidth',3);
hold on;
set(gca,'fontsize',20);
% title('Force Profile', 'Interpreter','LaTex','FontSize',20)
xlabel('$t [s]$','Interpreter','LaTex','FontSize',24);
ylabel('$\parallel F \parallel [N]$','Interpreter','LaTex','FontSize',24);
box on;
xlim([0,12])

figure;
plot(t, vecnorm(tau_ext_sedscrash_filt(1:6000,1:3)'),'linewidth',3);
hold on;
set(gca,'fontsize',20);
% title('Force Profile', 'Interpreter','LaTex','FontSize',24)
xlabel('$t [s]$','Interpreter','LaTex','FontSize',24);
ylabel('$\parallel F \parallel [N]$','Interpreter','LaTex','FontSize',24);
box on;
xlim([0,12])

%% plot effective tube on trajectory
target = [0;0];
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
x_ = x-repmat(target,1,size(x,2));
for i = 1:size(x_,2)
    xd(:,i) = feval(myds, x_(:,i));
end
% xd = feval(myds, x_);
omega_max = zeros(1,nx*ny);
delta = x_len * sigmascale;
for i = 1:size(x,2)
    omega_max(i) = pos_check(x_(:,i), x_cen, delta, B_);
%     if omega_max(i) > 2
%         omega_max(i) = 0;
%     end
end
for i = 1:size(x_,2)
    if omega_max(i) < 2
        xd(:,i) = feval(global_ds, x_(:,i));
    end
end
z_tmp = reshape(omega_max,nx,ny);
hcolor = pcolor(x_tmp,y_tmp,reshape(omega_max,nx,ny));
set(hcolor,'linestyle','none');
% load whiteCopperColorMap;
% colormap(flipud(cm));
% colormap(cubehelix(128,2.3,-0.61,1.97,0.75,[0.50,0.96],[0.20,0.90]));
colormap(cubehelix([],1.53,-1.32,2.74,0.81,[0.27,1],[0.38,0.86])); % temporary good parameter
% cubehelix_view(h_act)
colorbar;
caxis([min(min(z_tmp)), max(max(z_tmp))]);
h = streamslice(x_tmp,y_tmp,reshape(xd(1,:),ny,nx),reshape(xd(2,:),ny,nx),1,'method','cubic');
set(h,'LineWidth', 1)
set(h,'color',[0.0667  0.0667 0.0667]);
xl = xlabel('$x_y [m]$');
yl = ylabel('$x_z [m]$ ');
set(gca,'fontsize',20,'LineWidth',1);
set([xl yl],'interpreter','Latex','fontsize',24);
xlim([axlim(1),axlim(2)])
ylim([axlim(3),axlim(4)])
box on;
% set(gcf, 'Renderer', 'opengl')
axis(limits);
scatter(att(1),att(2), 150, [0 0 0],'d','Linewidth',2); 

%% plot sampled points
h_act = figure(); hold on;
length_seq = size(x_rec,2);
plot(x_rec(1,1),x_rec(2,1),'b*','markersize',10,'Linewidth',2);
plot(x_rec(1,2:length_seq),x_rec(2,2:length_seq),'r.','markersize',33);
box on;
scatter(att(1),att(2), 200, [0 0 0],'d','Linewidth',2); hold on;
axis(limits)
[hds_rob1] = plot_ds_model(h_act, global_ds, [0;0], limits,'medium'); hold on;
set(gca,'Fontsize',20);
xlabel('$x_y [m]$','Interpreter','LaTex','FontSize',24);
ylabel('$x_z [m]$','Interpreter','LaTex','FontSize',24);
% set(gca,'XTick',[], 'YTick', [])
plot(x_rec(1,1),x_rec(2,1),'b*','markersize',10,'Linewidth',2);
plot(x_rec(1,2:length_seq),x_rec(2,2:length_seq),'r.','markersize',33);
lgd=legend({'Starting Point', 'Sequence Points','Goal Point'},'Interpreter','LaTex','FontSize',16)
% lgd.EdgeColor=colorchange('d3d3d3');
% lgd.Color=colorchange('d3d3d3');
% title('Intepolated Data Points', 'Interpreter','LaTex','FontSize',20)
hold on;
%% plot stiffness in insertion task
t = -0.02:0.001:0.179;
m = (t+0.02)/(0.179+0.02);
k22 = zeros(1,length(t));
k11 = zeros(1,length(t));
t_max = 0.05;
t_min = 0.0;
k22_min = 100;
k22_diff = 600;
e_max = 0.05;
e_min = 0.0;
k11_min = 700;
k11_diff = 500;
for i = 1:length(t)
    k22(i) = k22_min + smooth_transition_rising(t(i),t_max,t_min,k22_diff);
    k11(i) = k11_min + k11_diff*smooth_transition_fall(t(i),e_max,e_min);
end

figure;
hold on;
plot(m,flip(k11),'r','LineWidth',4);
plot(m,flip(k22),'g','LineWidth',4);

set(gca,'Fontsize',20);
box on
xlabel('$s$','FontSize',24,'Interpreter','latex')
ylabel('$K \; [N/m]$','FontSize',24,'Interpreter','latex')
xlim([0 1])
ylim([0 1300])
legend({'$k_{1}$','$k_{2}$'},'Interpreter','latex','FontSize',20)

%% plot insertion task
X_rob_1=load('x_rob_socket1.txt') ;
X_rob_2=load('x_rob_socket2.txt') ;
X_rob_3=load('x_rob_socket3.txt') ;
X_rob_4=load('x_rob_socket4.txt') ;

fs = 500;
[B,A]=butter(1,1.2/(fs/2));

X_rob_1=filtfilt(B,A,X_rob_1);
X_rob_2=filtfilt(B,A,X_rob_2);
X_rob_3=filtfilt(B,A,X_rob_3);
X_rob_4=filtfilt(B,A,X_rob_4);

figure;
hold on;
plot(X_rob_1(:,2),X_rob_1(:,3),'-r','LineWidth',4);
plot(X_rob_2(:,2),X_rob_2(:,3),'-g','LineWidth',4);
plot(X_rob_3(:,2),X_rob_3(:,3),'-.b','LineWidth',4);
plot(X_rob_4(:,2),X_rob_4(:,3),':m','LineWidth',4);

set(gca,'Fontsize',20);
box on
xlabel('$x_y [m]$','FontSize',24,'Interpreter','latex')
ylabel('$x_z [m]$','FontSize',24,'Interpreter','latex')
xlim([0.472 0.502])
% legend({'$motion_{without,1}$','$motion_{without,2}$','$motion_{with,1}$','$motion_{with,2}$'},'Interpreter','latex')
legend({sprintf('motion 1\nno pert.'),sprintf('motion 2\nno pert.'),sprintf('motion 3\nwith pert.'),sprintf('motion 4\nwith pert.')},'FontSize',20,'Interpreter','latex')

%% plot motion with perturbation and motion regeneration
X_rob=load('x_rob_reg2.txt');  % suffix "_per1", "_SEDSper", "_reg1", "_reg2"
% limits = limits - [0 -0.05 -0.1 0];
fs = 500;
[B,A]=butter(1,1.2/(fs/2));
X_rob=filtfilt(B,A,X_rob);

figure;
hold on;
plot(x_rec(1,1:end), x_rec(2,1:end), 'r-.', 'linewidth', 4);
plot(X_rob(:,2),X_rob(:,3),'b','LineWidth',4);
scatter(att(1),att(2), 150, [0 0 0],'d','Linewidth',2); 

target = [0;0];
quality='high';
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
x_ = x-repmat(target,1,size(x,2));
xd = feval(global_ds, x_);
h = streamslice(x_tmp,y_tmp,reshape(xd(1,:),ny,nx),reshape(xd(2,:),ny,nx),1,'method','cubic');
set(h,'LineWidth', 1)
set(h,'color',[0.0667  0.0667 0.0667]);
xl = xlabel('$x_y [m]$ ');
yl = ylabel('$x_z [m]$ ');
set(gca,'fontsize',20);
set([xl yl],'interpreter','Latex','fontsize',24);
xlim([axlim(1),axlim(2)])
ylim([axlim(3),axlim(4)])
box on;
axis(limits);
legend({'reference','robot motion'},'Interpreter','latex')

%%
function result = smooth_transition_rising(t, t_max, t_min, scale)

	alpha_min = 0.0;

    if t>t_max
		alpha = 0;
    elseif t<t_min
		alpha = 1;
	else
		alpha = (0.5*(1 + cos(pi / (t_max - t_min)*(t - t_min)))*(1 - alpha_min) + alpha_min);
    end
    
    alpha = 1 - alpha;
    result = scale*alpha;

end

%%
function result = smooth_transition_fall(E, E_max, E_min)
	alpha_min = 0.0;

	if (E>E_max)
		alpha = 0;
    elseif (E<E_min)
		alpha = 1;
	else
		alpha = (0.5*(1 + cos(pi / (E_max - E_min)*(E - E_min)))*(1 - alpha_min) + alpha_min);
    end

    result=  alpha;
end
%% color change function
function color = colorchange(colorhex)
    color1 = hex2dec(colorhex(1:2))/255;
    color2 = hex2dec(colorhex(3:4))/255;
    color3 = hex2dec(colorhex(5:6))/255;
    color = [color1 color2 color3];
end