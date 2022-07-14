%% org_streamlines
h_act = figure();
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
att = [0.4; 0.1];
ref_traj = load('org_success_traj_1/vsds.txt');
plot(ref_traj(:,1), ref_traj(:,2), 'r', 'LineWidth', 4, 'LineStyle','--');
hold on
plot(border_1(:,1),border_1(:,2),'m', 'LineWidth',4, 'LineStyle', ':');
hold on
scatter(att(1),att(2), 150, [0 0 0],'d','Linewidth',2); hold on;
plot(border_2(:,1),border_2(:,2),'m', 'LineWidth',4, 'LineStyle', ':');
hold on
 
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
xlabel('$x_y [m]$','Interpreter','LaTex','FontSize',20);
ylabel('$x_z [m]$','Interpreter','LaTex','FontSize',20);
legend({'$reference$', '$border$', '$target$'}, 'Interpreter','LaTex', 'FontSize',20);
% title(['Streamlines of GP-based LMDS'], 'Interpreter','latex', 'FontSize',20);

box on;
% legend({'$robot\,motion$', '$reference$', '$target$', '$wall$', '$border$'}, 'Interpreter','LaTex', 'FontSize',12);
axis(limits);

%% streamlines of refinement
new_pos = load('new_points_case1/pos_train.txt');
new_vel = load('new_points_case1/vel_train.txt');
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
ds.limits = [-0.2, 0.7, -0.2, 0.7];

A  = -0.4 * eye(2);
% myds = @(x) linear_dynamics(x, A_hat, att);  % linear dynamics
myds =  @(x) gp_ds (x, input, att);

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

ref_traj = load('success_traj_case1/vsds.txt');
plot(ref_traj(:,1), ref_traj(:,2), 'r', 'LineWidth', 4, 'LineStyle','--');
hold on;

plot(border_1(:,1),border_1(:,2),'m', 'LineWidth',4, 'LineStyle', ':');
hold on
scatter(att(1),att(2), 150, [0 0 0],'d','Linewidth',2); hold on;
plot(border_2(:,1),border_2(:,2),'m', 'LineWidth',4, 'LineStyle', ':');
hold on

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
xlabel('$x_y [m]$','Interpreter','LaTex','FontSize',20);
ylabel('$x_z [m]$','Interpreter','LaTex','FontSize',20);
legend({'$reference$', '$border$', '$target$'}, 'Interpreter','LaTex', 'FontSize',20);
% title(['Streamlines of GP-based LMDS'], 'Interpreter','latex', 'FontSize',20);

box on;
axis(limits);


%% streamlines of adaption
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

ds.limits = [-0.2, 0.7, -0.2, 0.7];
% ds.limits = [-0.3, 1.2, -0.4, 1.1]; % from -1 to 1 in both two dimensions
% att = [0; 0];  % the global attractor is at (0,0).
A  = -0.4 * eye(2);
% myds = @(x) linear_dynamics(x, A_hat, att);  % linear dynamics
myds =  @(x) gp_ds (x, input, att);

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

ref_traj = load('success_traj_case2/vsds.txt');
plot(ref_traj(:,1), ref_traj(:,2), 'r', 'LineWidth', 4, 'LineStyle','--');
hold on;
plot(border_1(:,1),border_1(:,2),'m', 'LineWidth',4, 'LineStyle', ':');
hold on
scatter(att(1),att(2), 150, [0 0 0],'d','Linewidth',2); hold on;
plot(border_2(:,1),border_2(:,2),'m', 'LineWidth',4, 'LineStyle', ':');
hold on

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
xlabel('$x_y [m]$','Interpreter','LaTex','FontSize',20);
ylabel('$x_z [m]$','Interpreter','LaTex','FontSize',20);
% title(['Streamlines of GP-based LMDS'], 'Interpreter','latex', 'FontSize',20);

box on;
axis(ds.limits);
legend({'$reference$', '$border$', '$target$'}, 'Interpreter','LaTex', 'FontSize',20);

%% plot after incremental learning for the case 'close to demonstration'.
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
selected = [];
for i = 1:160
    if mod(i,2) == 0
        selected(end+1,:) = new_pos(i,1:2);
    end
end
selected = [selected; new_pos(160:end,1:2)];
scatter(selected(:,1), selected(:,2), 'b'); hold on;

ds.limits = [-0.2, 0.7, -0.2, 0.7];
% ds.limits = [-0.3, 1.2, -0.4, 1.1]; % from -1 to 1 in both two dimensions
% att = [0; 0];  % the global attractor is at (0,0).
A  = -0.4 * eye(2);
% myds = @(x) linear_dynamics(x, A_hat, att);  % linear dynamics
myds =  @(x) gp_ds (x, input, att);
[hatt_rob] = scatter(att(1),att(2), 150, [0 0 0],'d','Linewidth',2); hold on;
% [0 0 0] is the color RGB setting. 'd' represents the type of the marker.
[hds_rob1] = plot_ds(h_act, myds, att, ds.limits,'medium'); hold on;
set(gca,'fontsize',20,'LineWidth',1);
xlabel('$x_y [m]$','Interpreter','LaTex','FontSize',20);
ylabel('$x_z [m]$','Interpreter','LaTex','FontSize',20);
hold on
box on
axis(ds.limits);
legend({'$selected\, points$', '$target$'}, 'Interpreter','Latex', 'FontSize',20);
%  scatter(new_pos(:,1), new_pos(:,2), 'm');

%% plot after incremental learning for the case 'far away from demonstration'.
% after demonstration, plot the new streamlines
new_pos = load('new_points_case1/pos_train.txt');
new_vel = load('new_points_case1/vel_train.txt');
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
selected = [];
for i = 121:201
    if mod(i,2) == 0
        selected(end+1,:) = new_pos(i,1:2);
    end
end
selected = [selected; new_pos(1:120,1:2)];
scatter(selected(:,1), selected(:,2), 'b'); hold on;

ds.limits = [-0.2, 0.7, -0.2, 0.7];
% ds.limits = [-0.3, 1.2, -0.4, 1.1]; % from -1 to 1 in both two dimensions
% att = [0; 0];  % the global attractor is at (0,0).
A  = -0.4 * eye(2);
% myds = @(x) linear_dynamics(x, A_hat, att);  % linear dynamics
myds =  @(x) gp_ds (x, input, att);
[hatt_rob] = scatter(att(1),att(2), 150, [0 0 0],'d','Linewidth',2); hold on;
% [0 0 0] is the color RGB setting. 'd' represents the type of the marker.
[hds_rob1] = plot_ds(h_act, myds, att, ds.limits,'medium'); hold on;
set(gca,'fontsize',20,'LineWidth',1);
xlabel('$x_y [m]$','Interpreter','LaTex','FontSize',20);
ylabel('$x_z [m]$','Interpreter','LaTex','FontSize',20);
hold on
box on
axis(ds.limits);
legend({'$selected\, points$', '$target$'}, 'Interpreter','Latex', 'FontSize',20);



