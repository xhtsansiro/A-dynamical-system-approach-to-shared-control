function Stats_Forceder=bar_plot(Data,var_name,conditions) ;


    figure
    hold on
    variab=Data;
    mean_vel_ac_sub=mean(variab,1) ;
    std_vel=std(variab,0,1) ;
    [p,tbl,stats] = anova1(variab,[],'off');
    % Stats_Forceder=multcompare(stats,'Display','off')  ;
    Stats_Forceder=RepeatedAnova(Data,var_name,conditions) ;

%     title(strcat(var_name, ' p=',num2str(p)),'interpreter','latex','fontsize',18,'fontweight','normal') ;
%     title(strcat(var_name),'interpreter','latex','fontsize',18,'fontweight','normal') ;

    box on;grid on;hold on;
    set(gca,'xticklabel',conditions)
    set(gca,'TickLabelInterpreter','latex','fontsize',12)
    bar(mean_vel_ac_sub)
    er=errorbar(1:length(conditions),mean_vel_ac_sub,std_vel) ;er.LineStyle = 'none';  er.Color = [0 0 0];
    er.LineWidth=1.5 ;xticks(1:length(conditions)) ;


end
