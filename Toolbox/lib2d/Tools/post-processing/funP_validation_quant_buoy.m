function funP_validation_quant_buoy(simuldata,measdata,buoy,setting,h)
output=simuldata.output;
Proj=simuldata.Proj;
X=output.X;Y=output.Y;
T=output.time;
indX=funC_closest(X,buoy.x);
indY=funC_closest(Y,buoy.y);
indT1=funC_closest(T,buoy.tinterv(1));
indT2=funC_closest(T,buoy.tinterv(2));

Measdata=measdata(3:end,2:end);
Tmeas=measdata(3:end,1);

  
    for ii=1:length(indX)
        EtaS=squeeze(output.eta(:,indY(ii),indX(ii)));
        EtaM=interp1(Tmeas,Measdata(:,ii),T,'spline');
        if setting.cb.timeshift==1
                if setting.cb.timeshift_def==1
                    [Xcorel,lags]=xcorr(EtaS-mean(EtaS),EtaM-mean(EtaM),'Coeff');
                    MaxCor=max(Xcorel);
                    indMaxCorel=closest(Xcorel,MaxCor);
                    ShiftOpt=lags(indMaxCorel);
                else
                   ShiftOpt=setting.cb.timeshift_dt_edit; 
                end
            EtaM = circshift(EtaM,ShiftOpt);
        end
        EtaS=EtaS(indT1:indT2);
        EtaM=EtaM(indT1:indT2);
        EtaS=EtaS-mean(EtaS);
        EtaM=EtaM-mean(EtaM);
        
        varEtaS=var(EtaS);
        varEtaM=var(EtaM);
        %Corr=(EtaS'*EtaM)/sqrt(varEtaM*varEtaS)/(length(EtaM));
        Corr=(EtaS'*EtaM)/norm(EtaS)/norm(EtaM);
        
        HsS=4*sqrt(varEtaS)*1000;
        HsM=4*sqrt(varEtaM)*1000;        
        MaxEtaS=max(EtaS)*1000;
        MaxEtaM=max(EtaM)*1000;
        MinEtaS=min(EtaS)*1000;
        MinEtaM=min(EtaM)*1000;
        MeanEtaS=mean(EtaS)*1000;
        MeanEtaM=mean(EtaM)*1000;
        TpS=funSP_find_peak_period(T(indT1:indT2),EtaS)*1000;
        TpM=funSP_find_peak_period(T(indT1:indT2),EtaM)*1000;
        Tm01S=funSP_find_mean_period_Tm01(T(indT1:indT2),EtaS)*1000;
        Tm01M=funSP_find_mean_period_Tm01(T(indT1:indT2),EtaM)*1000;
        SkS=funSP_find_skewness(EtaS)*1000;
        SkM=funSP_find_skewness(EtaM)*1000;
        AsS=-funSP_find_asymmetry(EtaS)*1000;
        AsM=-funSP_find_asymmetry(EtaM)*1000;
        KurtS=kurtosis(EtaS)*1000;
        KurtM=kurtosis(EtaM)*1000;
        
        tabdata(ii,:)={Corr,HsS/HsM,MaxEtaS/MaxEtaM,MinEtaS/MinEtaM,...
            MeanEtaS/MeanEtaM,TpS/TpM,Tm01S/Tm01M,SkS/SkM,AsS/AsM,KurtS/KurtM};
    end
     set(h.pp_validation_quant_buoy_table,'data',tabdata); 
    
if buoy.savedata==1
      quantdata=cell2mat(tabdata);
      save([Proj.workdir,'PP_validation_quant_buoy_data.txt'],'quantdata','-ascii');   
end
end