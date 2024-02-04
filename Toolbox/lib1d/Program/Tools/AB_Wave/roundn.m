function var=roundn(varin,n)
%% round to n decimal (n<0)
var=round(varin*10^(-n))*10^(n);

