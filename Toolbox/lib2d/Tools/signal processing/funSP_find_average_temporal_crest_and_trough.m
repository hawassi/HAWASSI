function [ATC,ATT]=funSP_find_average_temporal_crest_and_trough(T,var_t)
mean_var=mean(var_t);
datvar=[T var_t];

try
    [~, ind] = dat2wa(datvar,mean_var,'c2c');
    if ~isempty(ind)
        if ind>0
            var_t_TC=var_t(ind);
            var_t_C=var_t_TC(var_t_TC>mean_var);
            var_t_T=var_t_TC(var_t_TC<mean_var);
            ATC=mean(var_t_C);
            ATT=mean(var_t_T);
        else
            ATC=0;
            ATT=0;
        end
    else
        ATC=0;
        ATT=0;
    end
catch
    ATC=0;
    ATT=0;
end
end