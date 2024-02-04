function funP_validation_plot_buoy(spatsnap,simuldata,measdata,setting,axesfig,hfig,IndPlot)
varId=setting.var;
output=simuldata.output;
Proj=simuldata.Proj;
X=output.X;Y=output.Y;
T=output.time;
if setting.coarse.check==0
    stepT=1;
else
    stepT=setting.coarse.param;
end
CheckId.SaveFig=setting.savefig.check;
CheckId.SaveFig_type=setting.savefig.format;


indX=funC_closest(X,spatsnap.x_snap(IndPlot));
indY=funC_closest(Y,spatsnap.y_snap(IndPlot));
indtS1=funC_closest(T,spatsnap.tinterv(1));
indtS2=funC_closest(T,spatsnap.tinterv(2));
Tmeas=measdata(3:end,1);
MeasData=measdata(3:end,IndPlot+1);
var_xyM=interp1(Tmeas,MeasData,T,'spline');


cla(axesfig);

        
        if varId==1
            var_xyS=squeeze(output.eta(:,indY,indX));
            varnameY='elevation [m]';
            savename='Elevation';
        elseif varId==2
            if simuldata.model.phiForm==1
            var_xyS=squeeze(output.phi(:,indY,indX));
            varnameY='potential';
            savename='Potential';
            else
            var_xyS=squeeze(output.u(:,indY,indX));
            varnameY='velocity (in x dir.)';
            savename='VelocityInX';    
            end
        else
            var_xyS=squeeze(output.v(:,indY,indX));
            varnameY='velocity (in y dir.)';
            savename='VelocityInY';  
        end
           
            if setting.cb.timeshift==1
                if setting.cb.timeshift_def==1
                    [Xcorel,lags]=xcorr(var_xyS-mean(var_xyS),var_xyM-mean(var_xyM),'Coeff');
                    MaxCor=max(Xcorel);
                    indMaxCorel=closest(Xcorel,MaxCor);
                    ShiftOpt=lags(indMaxCorel);
                else
                   ShiftOpt=setting.cb.timeshift_dt_edit; 
                end
            var_xyM = circshift(var_xyM,ShiftOpt);
            end
           var_xyM=var_xyM(indtS1:stepT:indtS2);
           var_xyS=var_xyS(indtS1:stepT:indtS2);
            
        if setting.cb.signal==1
            horzAx=T(indtS1:stepT:indtS2);
            xlab='time [s]';
            titletext='';
            horznmin=horzAx(1);horznmax=horzAx(end);
        else
            omsig=funC_freqspace(T(indtS1:stepT:indtS2));
           
            varSigS=var(var_xyS);
            S_hatS=abs(fft(var_xyS)).^2;
           if setting.spsmooth.check==1
               S_hatS=smooth(S_hatS,setting.spsmooth.param); 
            end
            var0S=trapz(omsig(1:floor(end/2)),S_hatS(1:floor(end/2)));
            S_hatS=S_hatS*varSigS/var0S;
            
            varSigM=var(var_xyM);
            S_hatM=abs(fft(var_xyM)).^2;
             if setting.spsmooth.check==1
               S_hatM=smooth(S_hatM,setting.spsmooth.param); 
            end
            var0M=trapz(omsig(1:floor(end/2)),S_hatM(1:floor(end/2)));
            S_hatM=S_hatM*varSigM/var0M;
            
            
            if  setting.cb.spectrum_var==1
                varnameY='E(\omega) [m^2 s/rad]';
                xlab='\omega [rad/s]';
                var_xyS=S_hatS(1:floor(end/2));
                var_xyM=S_hatM(1:floor(end/2));
                horzAx=omsig(1:floor(end/2));
            else
                varnameY='E(f) [m^2/hz]';
                xlab='f [hz]';
                var_xyS=S_hatS(1:floor(end/2)).*(2*pi);
                var_xyM=S_hatM(1:floor(end/2)).*(2*pi);
                horzAx=omsig(1:floor(end/2))./(2*pi);
            end
            savename='VarianceDensitySpectrum';
            titletext='Variance density spectrum';
           
            horznmin=omsig(1);horznmax=omsig(floor(end/2));
        end
        
        plot(axesfig,horzAx,var_xyM,'-b',horzAx,var_xyS,'--r');
        legend(axesfig,'Meas','Simul')
        title(titletext);
        
       
   
    if setting.horznlim.check==1
        horznmin=setting.horznlim.param(1);
        horznmax=setting.horznlim.param(2);
    end
    
    if setting.zlim.check==0
        maxvar=max(max(var_xyS),max(var_xyM));
        minvar=min(min(var_xyS),min(var_xyM));
        zlimm=[minvar maxvar];
    else
        zlimm=setting.zlim.param';
    end
     ylim(axesfig,zlimm);
     xlim(axesfig,[horznmin horznmax]);

xlabel(axesfig,xlab);ylabel(axesfig,varnameY);


axes_properties(axesfig,1);

if setting.savedat==1
    varr=[horzAx var_xyM var_xyS];
    if setting.cb.signal==1
    save([Proj.workdir,'PP_validation_Buoy_data_signal',savename,'.txt'],'varr','-ascii');    
    else
    save([Proj.workdir,'PP_validation_Buoy_data_spectrum',savename,'.txt'],'varr','-ascii');
    end
end

if CheckId.SaveFig==1
    if strcmp(CheckId.SaveFig_type,'.eps')
        ax_old = axesfig;
        f_new = figure('unit','normalized','position',[1.1 1.1 0.7 0.7]);
        ax_new = copyobj(ax_old,f_new);
        axes_properties(ax_new,1);
        set(f_new,'unit','normalized','position',[1.1 1.1 0.7 0.7])
        set(ax_new,'outerposition',[0 0 1 1]);
        legend(ax_new,'Meas','Simul')
        saveas(ax_new,cell2mat([Proj.workdir,'PP_validation_Buoy_',num2str(IndPlot),'_',savename,CheckId.SaveFig_type]),'epsc');
        close
    else
        ax_old = axesfig;
        f_new = figure('unit','normalized','position',[1.1 1.1 0.7 0.7]);
        ax_new = copyobj(ax_old,f_new);
        axes_properties(ax_new,1);
        set(f_new,'unit','normalized','position',[1.1 1.1 0.7 0.7])
        set(ax_new,'outerposition',[0 0 1 1]);
        legend(ax_new,'Meas','Simul')
        saveas(ax_new,cell2mat([Proj.workdir,'PP_validation_Buoy_',num2str(IndPlot),'_',savename,CheckId.SaveFig_type]))
        close;
    end
end

end