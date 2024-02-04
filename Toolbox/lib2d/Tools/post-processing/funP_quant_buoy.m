function funP_quant_buoy(simuldata,buoy,h)
output=simuldata.output;
Proj=simuldata.Proj;
X=output.X;Y=output.Y;
T=output.time;
indX=funC_closest(X,buoy.x);
indY=funC_closest(Y,buoy.y);
indT1=funC_closest(T,buoy.tinterv(1));
indT2=funC_closest(T,buoy.tinterv(2));

if buoy.combine==1
    tt=T(indT1:indT2);dt=T(2)-T(1);
    Nt=length(tt);
    Ndat=length(indX);
    Tend=tt(1)+(Ndat*Nt-1)*dt;
    TT=linspace(tt(1),Tend,Ndat*Nt);
    iter=1;
    Eta=zeros(length(TT),1);
    for ii=1:Ndat
        Eta(iter:iter+Nt-1)=squeeze(output.eta(indT1:indT2,indY(ii),indX(ii)));
        iter=iter+Nt;
    end
    
        Hs=4*sqrt(var(Eta));
        MaxEta=max(Eta);
        MinEta=min(Eta);
        MeanEta=mean(Eta);
        Tp=funSP_find_peak_period(TT,Eta);
        Tm01=funSP_find_mean_period_Tm01(TT,Eta);
        Sk=funSP_find_skewness(Eta);
        As=-funSP_find_asymmetry(Eta);
        Kurt=kurtosis(Eta);
    tabdata={Hs, MaxEta, MinEta, MeanEta,Tp,Tm01, Sk,As,Kurt};
    set(h.pp_quant_buoy_table,'data',tabdata);       
else
  
    for ii=1:length(indX)
        Eta=squeeze(output.eta(indT1:indT2,indY(ii),indX(ii)));
        Hs=4*sqrt(var(Eta));
        MaxEta=max(Eta);
        MinEta=min(Eta);
        MeanEta=mean(Eta);
        Tp=funSP_find_peak_period(T,Eta);
        Tm01=funSP_find_mean_period_Tm01(T,Eta);
        Sk=funSP_find_skewness(Eta);
        As=-funSP_find_asymmetry(Eta);
        Kurt=kurtosis(Eta-mean(Eta));
        tabdata(ii,:)={Hs,MaxEta,MinEta,MeanEta,Tp,Tm01,Sk,As,Kurt};
    end
     set(h.pp_quant_buoy_table,'data',tabdata); 
end    
if buoy.savedata==1
      quantdata=cell2mat(tabdata);
      save([Proj.workdir,'PP_quant_buoy_data.txt'],'quantdata','-ascii');   
end
end