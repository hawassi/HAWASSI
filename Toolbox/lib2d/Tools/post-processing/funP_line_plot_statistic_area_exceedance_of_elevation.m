function funP_line_plot_statistic_area_exceedance_of_elevation(simuldata,spatsnap,setting,axesfig,h)
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

xlab='\eta/H_{s}ref [-]';
varPnow=varP(indt1:indt2,indy1:indy2,indx1:indx2);
etaref=linspace(min(min(min(varPnow))),max(max(max(varPnow))),200);


dts=T(2)-T(1);
Ntt=floor((T(indt2)-T(indt1))/dts/stepT);

cla(axesfig);


[jProgressBar,statusbarObj]=funGui_JavProgressBar(h);
set(jProgressBar,'Maximum',Ntt, 'Value',0);
jProgressBar.setStringPainted( true );
statusbarObj.setText('Calculating...');
[jbStop]=funGui_Java_stopbutton(statusbarObj);
jProgressBar.setStringPainted( true );
tic;
iter=1;
Exceed=zeros(Ntt,length(etaref));

for ii=1:Ntt
    IDstop=funGui_eventLoopStop(jbStop);
    if IDstop==1, break;end;
    
    
indtti=indt1+stepT*(ii-1);
if indtti>=indt2, break; end;
varTi=squeeze(varP(indtti,indy1:indy2,indx1:indx2));
for jj=1:length(etaref)
        Exceed(ii,jj)=length(varTi(varTi>etaref(jj))) ;    
end

if mod(iter,floor(0.1*Ntt))==0 || iter==floor(0.01*Ntt)
        set(jProgressBar,'Maximum',Ntt, 'Value',iter);
        ETA=funGui_remain_time(iter,Ntt);
        statusbarObj.setText(['time remaining=', num2str(ETA)]);
end
    iter=iter+1;
end
Ntt=iter-1;
Ntotal=((indx2-indx1)+1)*((indy2-indy1)+1)*Ntt;

dx=X(2)-X(1);dy=Y(2)-Y(1);dt=T(2)-T(1);
Prob_Areaexceed=sum(Exceed)./Ntotal;
set(jProgressBar,'Maximum',Ntt, 'Value',Ntt); 
jProgressBar.setVisible(0);
jbStop.setVisible(0);
statusbarObj.setText('done.');

%plot3(axesfig,Xmax,Ymax,Tnow,'or');

semilogy(axesfig,etaref/Hs,Prob_Areaexceed,'or')   
xlabel(axesfig,xlab)
ylabel(axesfig,'Area exceedance')

axes_properties(axesfig,1);
% 
% xlim(axesfig,[etaref(1) etaref(end)]);

 if setting.savedat==1
        varr=[(etaref/Hs)' Prob_Areaexceed' ];
        save([Proj.workdir,'PP_line_statistic_area_exceedance','.txt'],'varr','-ascii');
 end
    
if setting.zlim.check==1
ylim(axesfig,[setting.zlim.param(1) setting.zlim.param(2)]);
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
        saveas(ax_new,cell2mat([Proj.workdir,'PP_line_statistic_area_exceedance','_',varname,CheckId.SaveFig_type]),'epsc');
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
        saveas(ax_new,cell2mat([Proj.workdir,'PP_line_statistic_area_exceedance','_',varname,CheckId.SaveFig_type]))
        close;
    end
end

end