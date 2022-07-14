%% First check 
exp.trial ={'01','02', '03', '04', '05', '06', '07', '08', '09', '10'};
exp.name = {'Jinsong', 'Tianshu', 'Kun', 'Chaitanya', 'Farbod', 'Yiheng', 'Philipp', 'Xun', 'Bowen', 'Havertz'};
exp.condition = {'no_controller', 'impedance', 'flow', 'vsds'};
exp.dataterms = {'velocity.txt', 'robot_real.txt', 'inertia.txt', 'control_force.txt'};
% always follow the order 01: free mode, 02: impedance controller, 03: flow
% controller, 04: vsds controller. the number represents how many trials
% succeed.
exp.record = [1,3,1,3; 1,3,2,3; 3,3,2,3; 2,0,3,3; 2,3,2,2; 2,1,1,3; 0,2,3,2; 2,3,3,3; 1,2,2,3; 2,3,3,3];
% detailed record of the trials
exp.detailed(:,:,1) = [0, 1, 0, 1; 0, 1, 1, 1; 1, 1, 0, 1; 1, 0, 1, 1; 0, 1, 1, 1; 1, 1, 0, 1; 0, 1, 1, 0; 0, 1, 1, 1; 1, 1, 0, 1; 1, 1, 1, 1];
exp.detailed(:,:,2) = [0, 1, 1, 1; 0, 1, 1, 1; 1, 1, 1, 1; 0, 0, 1, 1; 1, 1, 0, 1; 1, 0, 1, 1; 0, 1, 1, 1; 1, 1, 1, 1; 0, 0, 1, 1; 0, 1, 1, 1];
exp.detailed(:,:,3) = [1, 1, 0, 1; 1, 1, 0, 1; 1, 1, 1, 1; 1, 0, 1, 1; 1, 1, 1, 0; 0, 0, 0, 1; 0, 0, 1, 1; 1, 1, 1, 1; 0, 1, 1, 1; 1, 1, 1, 1];   
% load the execution time from txt file 
exp.time = load('exe_time.txt');

% calculate successful rate of each conditions. (total 4 conditions)
exp.successrate = sum(exp.record,1)/(3*length(exp.trial));

%% Plot the successful rate of four conditions
X = categorical({'FreeMode','OpenLoop','Flow','VSDS'});
X = reordercats(X, {'FreeMode','OpenLoop','Flow','VSDS'});
% X = categorical({'Condition_1','Condition_2','Condition_3','Condition_4'});
% X = reordercats(X, {'Condition_1', 'Condition_2', 'Condition_3', 'Condition_4'});
figure()
b = bar(X, exp.successrate, 0.75);
xtips = b.XEndPoints; ytips = b.YEndPoints; labels = string(b.YData);
for i = 1:length(b.YData)
    labels(i) = strcat(num2str(roundn(b.YData(i)*100, -2)), '%');
end
% title('Jerk of the robot motion','interpreter','latex','fontsize',18)
%  set(gca,'fontsize',40,'LineWidth',1);
set(gca,'TickLabelInterpreter','latex','fontsize',20)
text(xtips, ytips, labels, 'HorizontalAlignment','center',...
    'VerticalAlignment','bottom', 'FontSize',16);
% title("Successful rate of execution", 'Interpreter','latex', 'FontSize',18)
grid on
ylim([0 1.1]);

%% Plot the average execution time
% first calculate average time of each person
time = [];
% normalize terms, if the trial fails, take the highest time of all trials
for i = 1:4
    for j = 1:30
        if exp.time(i,j) == 0
            if max(exp.time(i,3*ceil(j/3)-2: 3*ceil(j/3))) == 0
                exp.time(i,j) = max (exp.time(i,:));
            else
                exp.time(i,j) = max(exp.time(i,3*ceil(j/3)-2: 3*ceil(j/3)));
            end
        end
    end
end

for k = 1:10
    time(:,k) = sum(exp.time(:,3*k-2: 3*k), 2)/ 3;
end

time = time';  % tranverse it

conditions={'FreeMode','OpenLoop','Flow','VSDS'} ;
bar_plot(time,'ExecutionTime',conditions) ;
set(gca,'fontsize',20,'LineWidth',1);
ylabel('$[s]$','Interpreter','LaTex','FontSize',25);
% ylabel('$x_z [m]$','Interpreter','LaTex','FontSize',30);
% title('Task execution time ','interpreter','latex','fontsize',18)

%%

exp.time_average(1) = sum(exp.time(1,:))/ sum(exp.record(:,1));
exp.time_average(2) = sum(exp.time(2,:))/ sum(exp.record(:,2));
exp.time_average(3) = sum(exp.time(3,:))/ sum(exp.record(:,3));
exp.time_average(4) = sum(exp.time(4,:))/ sum(exp.record(:,4));
X = categorical({'FreeMode','ImpedanceController','VelocityController','VSDSController'});
Y = reordercats(X, {'FreeMode', 'ImpedanceController', 'VelocityController', 'VSDSController'});
b = bar(X, exp.time_average, 0.75);

conditions={'FreeMode','ImpedanceController','VelocityController','VSDSController'} ;
bar_plot(jerk_mean,'Jerk',conditions) ;
bar_plot(f_ext_mean,'ExternalForce',conditions) ;

xtips = b.XEndPoints; ytips = b.YEndPoints; labels = string(b.YData);
for i = 1:length(b.YData)
    labels(i) = strcat(num2str(b.YData(i)),'s');
end
text(xtips, ytips, labels, 'HorizontalAlignment','center',...
    'VerticalAlignment','bottom');
ylim([0 19]);

%% NASA TLX load
tlx_load = load("NASA_TLX.txt");

conditions={'FreeMode','OpenLoop','Flow','VSDS'} ;
bar_plot(tlx_load,'TaskLoad',conditions) ;
% title('Task load index','interpreter','latex','fontsize',18)
set(gca,'fontsize',20,'LineWidth',1);
ylim([0 100])

%% The motion itself
fs=500 ; dt=1/fs ;
[B,A]=butter(1,3/(fs/2)); 
% calculate the motion characteristics
for i = 1:length(exp.trial)
    for j = 1:4
   
%         if exp.record(i,j) ~= 0
%             for k =1: exp.record(i,j)
        for k=1:3
                path = strcat(exp.trial{1,i},'/',exp.name{i}, '_', exp.condition{j},'_', num2str(k),'_');
                % compute the jerk of the robot motion
                path_trj = strcat(path, exp.dataterms{2});
                trj = load(path_trj);   % if successful traj., load all
                selected = ceil(0.8*length(trj));
                if exp.detailed(i,j,k) == 0
                    trj = trj(1:selected,:);    % if not, load only part of the traj.
                end

                % filter the traj with butterworth filter
                trj_fil = filtfilt(B,A,trj);
                % calculate velocity and filter it
                v = diff(trj_fil)/dt;
                v_fil = filtfilt(B,A,v); 
                % calculate accerleration and filter it
                a = diff(v_fil)/dt;
                a_fil = filtfilt(B,A,a); 
                % calculate jerk and filter it
                jerk = diff(a_fil)/dt;
                jerk_fil = filtfilt(B,A,jerk); 
                magni_jerk = vecnorm (jerk_fil(:,2:3),2,2);
                exp.jerk{i,j}(1,k) = rms(magni_jerk);

                % calculate the force, first of all, load the velocity
                % load the velocity of haptic device, filter it, calculate
                % accerleration and filter it
                path_v_hd = strcat(path, exp.dataterms{1});
                v_hd = load(path_v_hd);
                if exp.detailed(i,j,k) == 0
                    v_hd = v_hd(1:selected,:);    % if not success, load only part of the traj.
                end

                v_hd_fil = filtfilt(B,A,v_hd);
                a_hd = diff(v_hd_fil)/dt;
                a_hd_fil = filtfilt(B,A,a_hd);
                
                % load the control force u_c
                path_c_f = strcat(path, exp.dataterms{4});
                c_f = load(path_c_f);
                if exp.detailed(i,j,k) == 0
                    c_f = c_f(1:selected,:);    % if not success, load only part of the traj.
                end
                new = zeros(size(c_f,1), 1);
                c_f = [new, c_f];
                % load and convert
                % inertia matrix
                path_inertia = strcat(path, exp.dataterms{3});
                inertia = load(path_inertia);
                if exp.detailed(i,j,k) == 0
                    inertia = inertia(1:selected,:);    % if not success, load only part of the traj.
                end

                ext_f = []; % reset the force for each trial of each person
                for x = 1: size(inertia, 1)-1
                    inertia_matrix{1,x} = [inertia(x,1), inertia(x,2), inertia(x,3);
                                           inertia(x,4), inertia(x,5), inertia(x,6);
                                           inertia(x,7), inertia(x,8), inertia(x,9)];
                    ext_f(x,:) = (inertia_matrix{1,x}*a_hd_fil(x,:)')' - c_f(x,:);
                end
                % ext_f is the vector of (Fx, Fy, Fz) of all timestamps in
                % one trajectory. we calculate the crossproduct of F and v,
                % the velocity of haptic device is v_hd_fil;, ext_f has 1
                % row less than v_hd_fil
                 power = []; energy = 0;
                 for ii = 1:length(ext_f)
                     power (ii,1) = abs (dot(ext_f(ii,2:3), v_hd_fil(ii,2:3)));
                     energy = energy + power(ii,1)*0.002; 
                 end
                 
                 exp.power{i,j}(1,k) = mean(power);   % save the mean power of this trial
                 exp.energy{i,j}(1,k) =  energy; 

                % the force magnitude 
                 magni_ext_f = vecnorm (ext_f(:,2:3),2,2);
%                 force = rms(magni_ext_f);
                 force = mean(magni_ext_f);
                 exp.f_ext{i,j}(1,k) = force;            
            end
%         else
%             exp.jerk{i,j} = 0;  
% %             exp.f_ext{i,j}(1,1) = 0; exp.f_ext{i,j}(1,2) = 0;  % x and y
%             exp.f_ext{i,j}(1,k) = 0;
%         end

    end
end

%% Do the Anova test to jerk and force
jerk_mean = []; f_ext_mean = [];
for i =1:10
    for j=1:4
        jerk_mean(i,j) = sum(exp.jerk{i,j})/3;
        f_ext_mean(i,j) = sum(exp.f_ext{i,j})/3;
        power_mean(i,j) = sum(exp.power{i,j})/3;
        energy_mean(i,j) = sum(exp.energy{i,j})/3;
    end
end
f = f_ext_mean;
f(:,2) =[];
% conditions={'FreeMode','Flow','VSDS'} ;
conditions={'FreeMode','OpenLoop','Flow','VSDS'} ;
bar_plot(jerk_mean,'Jerk',conditions) ;
ylabel('$[m/s^3]$','Interpreter','LaTex','FontSize',20);
% ylabel('$Jerk [m/s^3]$','Interpreter','LaTex','FontSize',20);
% title('Jerk of the robot motion','interpreter','latex','fontsize',18)
set(gca,'fontsize',20,'LineWidth',1);
ylim([0,0.6]);

% bar_plot(f,'ExternalForce',conditions) ;
% bar_plot(power_mean,'Power',conditions) ;
% bar_plot(energy_mean,'Energy',conditions) ;

%% Do the analysis of the questionnaire
q1 = load("question_1.txt");
q2 = load("question_2.txt");
q3 = load("question_3.txt");
q4 = load("question_4.txt");

average = []; deviation = [];  % represents mean and standard deviation
% conditions={'FreeMode','ImpedanceController','VelocityController','VSDSController'} ;
conditions={'OpenLoop','Flow','VSDS'} ;
% calculate the mean and std of the data in q2, q3, q4
average(1,:) = mean(q2(:,2:end));
average(2,:) = mean(q3(:,2:end));
average(3,:) = mean(q4(:,2:end));
deviation(1,:) = std(q2(:,2:end));
deviation(2,:) = std(q3(:,2:end));
deviation(3,:) = std(q4(:,2:end));

% plot the barplots

% X = categorical({'Guidance useful?','Fight Guidance?','Feel in control?'});
% X = reordercats(X,{'Guidance useful?','Fight Guidance?','Feel in control?'});

X = categorical({'Question 1','Question 2','Question 3'});
X = reordercats(X,{'Question 1','Question 2','Question 3'});
bar(X,average);

% put the error plot
hold on;
er = errorbar([1 2 3], average(:,2), deviation(:,2), 'k', 'LineStyle','None'); 
er.LineWidth = 1.5;
hold on 
er = errorbar([1.225 2.225 3.225], average(:,3), deviation(:,3), 'k', 'LineStyle','None');
er.LineWidth = 1.5;
hold on 
er = errorbar([0.775 1.775 2.775], average(:,1), deviation(:,1), 'k', 'LineStyle','None');
er.LineWidth = 1.5;

box on
set(gca,'TickLabelInterpreter','latex','fontsize',20)
legend({'$OpenLoop$', '$Flow$',  '$VSDS$' }, 'Interpreter','LaTex', 'FontSize',12);
ylim([0,6])
grid on

% bar_plot(q1,'FeelGuidance',conditions) ;

bar_plot(q2(:,2:end),'GuidanceUseful',conditions) ;
bar_plot(q3(:,2:end),'FightGuidance',conditions) ;
bar_plot(q4(:,2:end),'InControl',conditions) ;

% plot the bars of three questions in one figure.
