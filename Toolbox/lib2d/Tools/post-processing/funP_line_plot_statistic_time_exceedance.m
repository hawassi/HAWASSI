function funP_line_plot_statistic_time_exceedance(simuldata,spatsnap,setting,axesfig,h)
output=simuldata.output;
Proj=simuldata.Proj;
input=simuldata.input;
X=output.X;Y=output.Y;
dom=simuldata.dom;
T=output.time;
influx=simuldata.influx;
if isfield(simuldata,'ivp')
ivp=simuldata.ivp;
end
if isfield(simuldata,'bdyassim')
bdyassim=simuldata.bdyassim;
end

varP=output.eta;
varname='elevation';


if setting.coarse.check==0
stepT=1;
else
stepT=setting.coarse.param;  
end
CheckId.SaveFig=setting.savefig.check;
CheckId.SaveFig_type=setting.savefig.format;

if spatsnap.area_check==1
    indx1=funC_closest(X,spatsnap.area_x_interv(1));
    indx2=funC_closest(X,spatsnap.area_x_interv(2));
    indy1=funC_closest(Y,spatsnap.area_y_interv(1));
    indy2=funC_closest(Y,spatsnap.area_y_interv(2));
else
    indy1=1;indy2=length(Y);
    indx1=1;indx2=length(X);
end

if setting.tlim.check==1
indt1=funC_closest(T,setting.tlim.param(1));
indt2=funC_closest(T,setting.tlim.param(2));
else
indt1=1;indt2=length(T);   
end


if setting.Hsref.check==0
if strcmpi(input.wave.option,'Yes')
    Hs=max(influx.par.Hs);
else
    if ~strcmpi(ivp.typename,'Zero')
        Hs=ivp.par.Hs;
    elseif input.bdyassim.option==1
        Hs=bdyassim.par.Hs;
    else
        Hs=0;
    end
end
else
Hs=setting.Hsref.param;        
end



cla(axesfig);




xlab='\eta/H_s [-]';

varPnow=varP(indt1:indt2,indy1:indy2,indx1:indx2);
etaref=linspace(min(min(min(varPnow))),max(max(max(varPnow))),200)

Nx=((indx2-indx1)+1);
Ny=((indy2-indy1)+1);
Ntotal=Nx*Ny;%%point over area

Ntt=length(T(indt1:indt2));
Ndat=length(etaref);
timeExceed=zeros(Ndat,1);
iter=1;

% [jProgressBar,statusbarObj]=funGui_JavProgressBar(h);
% set(jProgressBar,'Maximum',Ndat, 'Value',0);
% jProgressBar.setStringPainted( true );
% statusbarObj.setText('Calculating...');
% [jbStop]=funGui_Java_stopbutton(statusbarObj);
% jProgressBar.setStringPainted( true );
% tic;


% for ii=1:Ndat
%     IDstop=funGui_eventLoopStop(jbStop);
%     if IDstop==1, break;end;
    
%    timeExceed(ii)=length(varPnow(varPnow>etaref(ii)))/Ntotal;
%    
% 
% if mod(iter,floor(0.1*Ndat))==0 || iter==floor(0.01*Ndat)
%         set(jProgressBar,'Maximum',Ndat, 'Value',iter);
%         ETA=funGui_remain_time(iter,Ndat);
%         statusbarObj.setText(['time remaining=', num2str(ETA)]);
% end
%     iter=iter+1;
% end
% Prob_TimeExceed=timeExceed./Ntt;

[jProgressBar,statusbarObj]=funGui_JavProgressBar(h);
set(jProgressBar,'Maximum',Ntotal, 'Value',0);
jProgressBar.setStringPainted( true );
statusbarObj.setText('Calculating...');
[jbStop]=funGui_Java_stopbutton(statusbarObj);
jProgressBar.setStringPainted( true );
tic;


timeExceed=zeros(Ndat,1);
for ii=1:Nx
    for jj=1:Ny
        IDstop=funGui_eventLoopStop(jbStop);
        if IDstop==1, break;end;
        indtti=indt1+stepT*(ii-1);
        if indtti>=indt2, break; end;
        varTi=squeeze(varP(indt1:indt2,indy1+jj-1,indx1+ii-1));
        
        for kk=1:Ndat
        timeExceed(kk)=timeExceed(kk)+length(varTi(varTi>etaref(kk)))./Ntt;
        end
        
        if mod(iter,floor(0.1*Ntotal))==0 || iter==floor(0.01*Ntotal)
            set(jProgressBar,'Maximum',Ntotal, 'Value',iter);
            ETA=funGui_remain_time(iter,Ntotal);
            statusbarObj.setText(['time remaining=', num2str(ETA)]);
        end
        iter=iter+1;
    end
end



Prob_TimeExceed=timeExceed./Ntotal;

set(jProgressBar,'Maximum',Ntotal, 'Value',Ntotal); 
jProgressBar.setVisible(0);
jbStop.setVisible(0);
statusbarObj.setText('done.');

%plot3(axesfig,Xmax,Ymax,Tnow,'or');

semilogy(axesfig,etaref/Hs,Prob_TimeExceed,'or')   
xlabel(axesfig,xlab)
ylabel(axesfig,'Time exceedance]')
 %  legend(axesfig,'+\eta','-\eta')
axes_properties(axesfig,1);
% 
% xlim(axesfig,[etaref(1) etaref(end)]);
 if setting.savedat==1
        varr=[(etaref/Hs)' Prob_TimeExceed'];
        save([Proj.workdir,'PP_line_statistic_time_exceedance','.txt'],'varr','-ascii');
 end
 
if setting.zlim.check==1
ylim(axesfig,[setting.zlim.param(1) setting.zlim.param(2)]);
else
% else
% ylim(axesfig,[min(min(Prob_TimeExceedmin),min(Prob_TimeExceedplus)) ...
%     max(max(Prob_TimeExceedmin),max(Prob_TimeExceedplus))+0.00001])
end

if setting.spatlim.check==1
xlim(axesfig,[setting.spatlim.param(1) setting.spatlim.param(2)]);
end
axes_properties(axesfig,1);


if CheckId.SaveFig==1
    if strcmp(CheckId.SaveFig_type,'.eps')
        ax_old = axesfig;
        f_new = figure('unit','normalized','position',[1.1 1.1 0.7 0.7]);
        ax_new = copyobj(ax_old,f_new);
        set(f_new,'visible','off')
        axes(ax_new)
        axes_properties(ax_new,1);
        set(f_new,'unit','normalized','position',[1.1 1.1 0.7 0.7])
        set(ax_new,'outerposition',[0 0 1 1]);
        saveas(ax_new,cell2mat([Proj.workdir,'PP_line_statistic_time_exceedance_',varname,CheckId.SaveFig_type]),'epsc');
        close
    else
        ax_old = axesfig;
        f_new = figure('unit','normalized','position',[1.1 1.1 0.7 0.7]);
        ax_new = copyobj(ax_old,f_new);
        set(f_new,'visible','off')
        axes(ax_new)
        axes_properties(ax_new,1);
        set(f_new,'unit','normalized','position',[1.1 1.1 0.7 0.7])
        set(ax_new,'outerposition',[0 0 1 1]);
        saveas(ax_new,cell2mat([Proj.workdir,'PP_line_statistic_time_exceedance_',varname,CheckId.SaveFig_type]))
        close;
    end
end

end