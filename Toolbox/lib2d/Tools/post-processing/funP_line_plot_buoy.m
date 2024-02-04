function funP_line_plot_buoy(spatsnap,simuldata,setting,axesfig,hfig)
cla;
varId=setting.var;
output=simuldata.output;
input=simuldata.input;
dom=simuldata.dom;
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

indX=funC_closest(X,spatsnap.x_snap);
indY=funC_closest(Y,spatsnap.y_snap);

if  setting.cb.spectrum_var<3
    if setting.cb.combine_buoys==1
        indtd1=funC_closest(T,setting.cb.combine_tinterv(1));
        indtd2=funC_closest(T,setting.cb.combine_tinterv(2));
        
        tt=output.time(indtd1:indtd2);dt=tt(2)-tt(1);
        Nt=length(tt);Ndat=length(indX);
        Tend=tt(1)+(Ndat*Nt-1)*dt;
        T=linspace(tt(1),Tend,Ndat*Nt);
        
        
        if setting.tlim.check==1
            indt1=funC_closest(T,setting.tlim.param(1));
            indt2=funC_closest(T,setting.tlim.param(2));
        else
            indt1=1;indt2=length(T);
        end
        
        var_xy=zeros(Nt*Ndat,1);
        iter=1;
        
        for ii=1:Ndat
            if varId==1
                var_xy(iter:iter+Nt-1)=squeeze(output.eta(indtd1:indtd2,indY(ii),indX(ii)));
            elseif varId==2
                if simuldata.model.phiForm==1
                var_xy(iter:iter+Nt-1)=squeeze(output.phi(indtd1:indtd2,indY(ii),indX(ii)));
                else
                var_xy(iter:iter+Nt-1)=squeeze(output.u(indtd1:indtd2,indY(ii),indX(ii)));
                end
            elseif varId==3
                var_xy(iter:iter+Nt-1)=squeeze(output.v(indtd1:indtd2,indY(ii),indX(ii)));
            end
            iter=iter+Nt;
        end
        if varId==1
            varnameY='elevation [m]';
            savename='Elevation';
        elseif varId==2
            if simuldata.model.phiForm==1
            varnameY='potential';
            savename='Potential';
            else
            varnameY='velocity (in x dir.)';
            savename='Velocity_inX';    
            end
        else
            varnameY='velocity (in y dir.)';
            savename='Velocity_inY'; 
        end
        xlab='time [s]';
        titletext='Combined signal';
        var_xy=var_xy(indt1:stepT:indt2);
        T=T(indt1:stepT:indt2);
        horzAx=T';
        
        
        if setting.cb.signal==0
            omsig=funC_freqspace(T);
            var_xy=var_xy-mean(var_xy);
            varSig=var(var_xy);
            S_hat=abs(fft(var_xy)).^2;
            if setting.spsmooth.check==1
                S_hat=smooth(S_hat,setting.spsmooth.param);
            end
            
            var0=trapz(omsig(1:floor(end/2)),S_hat(1:floor(end/2)));
            S_hat=S_hat*varSig/var0;
            
            
            if  setting.cb.spectrum_var==1
                varnameY='E(\omega) [m^2 s/rad]';
                xlab='\omega [rad/s]';
                var_xy=S_hat(1:floor(end/2));
                horzAx=omsig(1:floor(end/2));
            else
                varnameY='E(f) [m^2/hz]';
                xlab='f [hz]';
                var_xy=S_hat(1:floor(end/2)).*(2*pi);
                horzAx=omsig(1:floor(end/2))./(2*pi);
            end
            savename='VarianceDensitySpectrum';
            titletext='Variance density spectrum';
            
            
        end
        maxvar=max(var_xy);
        minvar=min(var_xy);
        
        plot(axesfig,horzAx,var_xy,'r');
        title(titletext);
        
    else
        if setting.tlim.check==1
            indt1=funC_closest(T,setting.tlim.param(1));
            indt2=funC_closest(T,setting.tlim.param(2));
        else
            indt1=1;indt2=length(T);
        end
        
        strc={'r','--g','-.b','--k','-.c'};
        indstr=1;
        hold(axesfig,'off');
        maxvar=0;minvar=0;
        for ii=1:length(indX)
            
            if setting.cb.signal==1
                if varId==1
                    var_xy=squeeze(output.eta(indt1:stepT:indt2,indY(ii),indX(ii)));
                    varnameY='elevation [m]';
                    savename='Elevation';
                elseif varId==2
                     if simuldata.model.phiForm==1
                        var_xy=squeeze(output.phi(indt1:stepT:indt2,indY(ii),indX(ii)));
                        varnameY='potential';
                        savename='Potential';
                    else
                        var_xy=squeeze(output.u(indt1:stepT:indt2,indY(ii),indX(ii)));
                        varnameY='velocity (in x dir.)';
                        savename='Velocity_inX';
                    end
                else
                        var_xy=squeeze(output.v(indt1:stepT:indt2,indY(ii),indX(ii)));
                        varnameY='velocity (in y dir.)';
                        savename='Velocity_inY';
                end
           
                
                horzAx=T(indt1:stepT:indt2);
                xlab='time [s]';
                titletext='';
            else
                Etaa=squeeze(output.eta(indt1:stepT:indt2,indY(ii),indX(ii)));
                omsig=funC_freqspace(T(indt1:stepT:indt2));
                Etaa=Etaa-mean(Etaa);
                varSig=var(Etaa);
                S_hat=abs(fft(Etaa)).^2;
                if setting.spsmooth.check==1
                    S_hat=smooth(S_hat,setting.spsmooth.param);
                end
                var0=trapz(omsig(1:floor(end/2)),S_hat(1:floor(end/2)));
                S_hat=S_hat*varSig/var0;
                
                
                if  setting.cb.spectrum_var==1
                    varnameY='E(\omega) [m^2 s/rad]';
                    xlab='\omega [rad/s]';
                    var_xy=S_hat(1:floor(end/2));
                    horzAx=omsig(1:floor(end/2));
                else
                    varnameY='E(f) [m^2/hz]';
                    xlab='f [hz]';
                    var_xy=S_hat(1:floor(end/2)).*(2*pi);
                    horzAx=omsig(1:floor(end/2))./(2*pi);
                end
                savename='VarianceDensitySpectrum';
                titletext='Variance density spectrum';
                
            end
            
            if indstr==6,indstr=1;end
            plot(axesfig,horzAx,var_xy,strc{indstr});
            title(titletext);
            hold(axesfig,'on');
            indstr=indstr+1;
            maxvar=max(maxvar,max(var_xy));
            minvar=min(minvar,min(var_xy));
        end
    end
    
    if setting.cb.signal==1
        if setting.zlim.check==0
            if setting.cb.combine_buoys==0
                if varId==1
                    maxvar=max(max(max(output.eta)));
                    minvar=min(min(min(output.eta)));
                elseif varId==2
                     if simuldata.model.phiForm==1
                    maxvar=max(max(max(output.phi)));
                    minvar=min(min(min(output.phi)));
                     else
                    maxvar=max(max(max(output.u)));
                    minvar=min(min(min(output.u)));
                     end
                else
                    maxvar=max(max(max(output.v)));
                    minvar=min(min(min(output.v)));
                end
            else
                maxvar=max(var_xy);
                minvar=min(var_xy);
            end
            
            
            zlimm=[minvar maxvar];
        else
            zlimm=setting.zlim.param';
        end
        ylim(axesfig,zlimm);
        if setting.horznlim.check==0
            xlimm=[T(indt1) T(indt2)];
        else
            xlimm=[setting.horznlim.param(1) setting.horznlim.param(2)];
        end
        
    else
        if setting.zlim.check==1
            zlimm=setting.zlim.param';
        else
            zlimm=[minvar maxvar];
        end
        ylim(axesfig,zlimm);
        
        if setting.horznlim.check==0
            xlimm=[horzAx(1) horzAx(end)];
        else
            xlimm=[setting.horznlim.param(1) setting.horznlim.param(2)];
        end
    end
    
    xlabel(axesfig,xlab);ylabel(axesfig,varnameY);
    
    xlim(axesfig,xlimm);
    axes_properties(axesfig,1);
    
else
    cla(axesfig,'reset');
    if setting.tlim.check==1
        indt1=funC_closest(T,setting.tlim.param(1));
        indt2=funC_closest(T,setting.tlim.param(2));
    else
        indt1=1;indt2=length(T);
    end
     pos=zeros(length(indX),3);
     meanDepth=0;
     timeS=T(indt1:stepT:indt2);
     datSignal=zeros(length(timeS),length(indX)+1);
     datSignal(:,1)=timeS;
     for ii=1:length(indX)
            
                if varId==1
                    var_xy=squeeze(output.eta(indt1:stepT:indt2,indY(ii),indX(ii)));
                    varnameY='elevation [m]';
                    savename='DirSpectElevation';
                elseif varId==2
                    if simuldata.model.phiForm==1
                    var_xy=squeeze(output.phi(indt1:stepT:indt2,indY(ii),indX(ii)));
                    varnameY='potential';
                    savename='DirSpectPotential';
                    else
                    var_xy=squeeze(output.u(indt1:stepT:indt2,indY(ii),indX(ii)));
                    varnameY='velocity (in x dir)';
                    savename='DirSpectVelocityInX';    
                    end
                else
                    var_xy=squeeze(output.v(indt1:stepT:indt2,indY(ii),indX(ii)));
                    varnameY='velocity (in y dir)';
                    savename='DirSpectVelocityInY';
                end
               pos(ii,1:2)=[X(indX(ii)) Y(indY(ii))];
               meanDepth=meanDepth+dom.bathy.profile(indY(ii),indX(ii));
               datSignal(:,ii+1)=var_xy;
     end            
    meanDepth=-meanDepth/length(indX);
    type = ones(1,length(indX)); bfs = ones(1,length(indX));
    
    S = dat2dspec(datSignal,[pos type' bfs'],meanDepth);
    var_xy=S.S;
    horzAx=S.w;
    if  setting.cb.spectrum_var==4
    S.f=S.w./2./pi;
    S=rmfield(S,'w');
    horzAx=S.f;
    end
    plotspec(S,axesfig);

end

if setting.savedat==1
    
    if setting.cb.signal==1
        varr=[horzAx var_xy];
        save([Proj.workdir,'PP_line_signal_data_',savename,'.txt'],'varr','-ascii');
    else
        varr=[horzAx var_xy];
        save([Proj.workdir,'PP_line_vardensspect_elevation_data.txt'],'varr','-ascii');
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
            saveas(ax_new,cell2mat([Proj.workdir,'PP_line_buoy_',savename,CheckId.SaveFig_type]),'epsc');
            close
        else
            ax_old = axesfig;
            f_new = figure('unit','normalized','position',[1.1 1.1 0.7 0.7]);
            ax_new = copyobj(ax_old,f_new);
            axes_properties(ax_new,1);
            set(f_new,'unit','normalized','position',[1.1 1.1 0.7 0.7])
            set(ax_new,'outerposition',[0 0 1 1]);
            saveas(ax_new,cell2mat([Proj.workdir,'PP_line_buoy',savename,CheckId.SaveFig_type]))
            close;
        end
    end

end