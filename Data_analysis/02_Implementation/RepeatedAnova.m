function tbl= RepeatedAnova(Data,Measure_name,conditions)


    % conditions={'Low','Var','EMG','High'} ;
    % Data=Jerk_Mot ;
    % 
    % Measure_name='JerkMotion' ;

    n_sub=length(Data) ;
    n_cond=length(conditions) ;
    counter_tot=1 ;

    for k=1:n_cond
        for i=1:length(Data)
            conditions_vec_tmp{counter_tot}=conditions{k} ;
            counter_tot=counter_tot+1 ;
        end
    end
    conditions_names_vec=conditions_vec_tmp' ;
    meas=Data(:) ;

    t = table(conditions_names_vec,meas(:,1),...
    'VariableNames',{'conditions',Measure_name});
    Meas = table([1]','VariableNames',{'Measurements'});

    resp_pred=strcat(Measure_name,'~conditions') ;

    rm = fitrm(t,resp_pred,'WithinDesign',Meas) ;


    [p,t,stats] =ranova(rm,'WithinModel',eye(1)) ;

    tbl = multcompare(rm,'conditions') ;

end