function corr=fun_corr(meas,simul)
corr=(simul'*meas)/sqrt(var(meas)*var(simul))/(length(meas));
if corr>1
    corr=1;
elseif corr<-1
    corr=-1;
end