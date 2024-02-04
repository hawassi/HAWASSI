function funP_line_plot_Extreme(simuldata,setting,axesfig,h)
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

if setting.xlim.check==1
indx1=funC_closest(X,setting.xlim.param(1));
indx2=funC_closest(X,setting.xlim.param(2));
else
indx1=1;indx2=length(X); 
end

if setting.ylim.check==1
indy1=funC_closest(Y,setting.ylim.param(1));
indy2=funC_closest(Y,setting.ylim.param(2));
else
indy1=1;indy2=length(Y); 
end  
if setting.tlim.check==1
indt1=funC_closest(T,setting.tlim.param(1));
indt2=funC_closest(T,setting.tlim.param(2));
else
indt1=1;indt2=length(T);   
end

if setting.ampliref.check==0
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
Hs=setting.ampliref.param;        
end
Xlabel='x[m]';
Ylabel='y[m]';
Zlabel='t[s]';
if setting.plotopt==1 
savevar='Extreme crests';
else
savevar='Extreme troughs';    
end

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
indx=1;
Xxtrem=[];Yxtrem=[];ExtVar=[];Tnow=[];
[XX,YY]=meshgrid(X(indx1:indx2),Y(indy1:indy2));

for ii=1:Ntt
    IDstop=funGui_eventLoopStop(jbStop);
    if IDstop==1, break;end;
    
indtti=indt1+stepT*(ii-1);
if indtti>=indt2, break; end;
varTi=squeeze(varP(indtti,indy1:indy2,indx1:indx2));
   
if setting.plotopt==1  
 Idx=find(varTi>Hs) ;
else
Idx=find(varTi<-Hs) ;
end
if ~isempty(Idx)
Xxtrem(indx:indx+length(Idx)-1)=XX(Idx);
Yxtrem(indx:indx+length(Idx)-1)=YY(Idx);
ExtVar(indx:indx+length(Idx)-1)=varTi(Idx);
Tnow(indx:indx+length(Idx)-1)=T(indtti);
indx=indx+length(Idx);
end

if mod(iter,floor(0.1*Ntt))==0 || iter==floor(0.01*Ntt)
        set(jProgressBar,'Maximum',Ntt, 'Value',iter);
        ETA=funGui_remain_time(iter,Ntt);
        statusbarObj.setText(['time remaining=', num2str(ETA)]);
end

    iter=iter+1;
end

set(jProgressBar,'Maximum',Ntt, 'Value',Ntt); 
jProgressBar.setVisible(0);
jbStop.setVisible(0);
statusbarObj.setText('done.');

%plot3(axesfig,Xmax,Ymax,Tnow,'or');


if ~isempty(ExtVar)
scatter3(axesfig, Xxtrem(:), Yxtrem(:), Tnow(:), [],ExtVar(:), 'filled' )
title(axesfig,savevar);
view([-40 10]);
colormap(axesfig,'jet');
axes_properties(axesfig,1);
axes(axesfig)
% cbfreeze;cblabel(varname);
% freezeColors;
cb=colorbar;
ylabel(cb,[varname, ' [m]'])
axes_properties_cb(cb,1.5);
zlabel(axesfig,Zlabel);
ylabel(axesfig,Ylabel);
xlabel(axesfig,Xlabel);
 set(h.monitorbox,'String',['>> '],'foregroundcolor','k');
else
 set(h.monitorbox,'String',['>>There is no extreme event.'],'foregroundcolor','k');
end

axes_properties(axesfig,1);

xlim(axesfig,[X(indx1) X(indx2)]);ylim(axesfig,[Y(indy1) Y(indy2)]);
zlim(axesfig,[T(indt1) T(indt2)]);

axes_properties(axesfig,1);


if CheckId.SaveFig==1
    if strcmp(CheckId.SaveFig_type,'.eps')
        ax_old = axesfig;
        f_new = figure('unit','normalized','position',[1.1 1.1 0.7 0.7]);
        ax_new = copyobj(ax_old,f_new);
        set(f_new,'visible','off')
        axes(ax_new)
        if setting.plotopt==4 || setting.plotopt==5
        colormap(ax_new,'jet');
        cb=colorbar;
        ylabel(cb,varname,'fontweight', 'bold')
        end
        axes_properties(ax_new,1);
        axes_properties_cb(cb,1);
        set(f_new,'unit','normalized','position',[1.1 1.1 0.7 0.7])
        set(ax_new,'outerposition',[0 0 1 1]);
        saveas(ax_new,cell2mat([Proj.workdir,'PP_line_',savevar,'_',varname,CheckId.SaveFig_type]),'epsc');
        close
    else
        ax_old = axesfig;
        f_new = figure('unit','normalized','position',[1.1 1.1 0.7 0.7]);
        ax_new = copyobj(ax_old,f_new);
        set(f_new,'visible','off')
        axes(ax_new)
        if setting.plotopt==4 || setting.plotopt==5
        colormap(ax_new,'jet');
        cb=colorbar;
        ylabel(cb,varname,'fontweight', 'bold')
        end
        axes_properties(ax_new,1);axes_properties(cb,1);
        set(f_new,'unit','normalized','position',[1.1 1.1 0.7 0.7])
        set(ax_new,'outerposition',[0 0 1 1]);
        saveas(ax_new,cell2mat([Proj.workdir,'PP_line_',savevar,'_',varname,CheckId.SaveFig_type]))
        close;
    end
end

end