function funP_density_plot_statistic(cbvar,simuldata,setting,axesfig,h)
hfig=h.fig;
varId=cbvar;
output=simuldata.output;
input=simuldata.input;
dom=simuldata.dom;
Proj=simuldata.Proj;
if isfield(simuldata,'bdyassim')
    bdyassim=simuldata.bdyassim;
end
xi=setting.xlim.param(1);
xf=setting.xlim.param(2);
yi=setting.ylim.param(1);
yf=setting.ylim.param(2);
ti=setting.tlim.param(1);
tf=setting.tlim.param(2);
if setting.view.check==1
    viewplot=setting.view.param';
else
    viewplot=2;
end
X=output.X;Y=output.Y;
T=output.time;
if setting.coarse.check==0
    stepS=1;stepT=1;
else
    stepT=setting.coarse.param;
    stepS=1;
end
CheckId.SaveFig=setting.savefig.check;
CheckId.SaveFig_type=setting.savefig.format;



if setting.xlim.check==1
    indx1=funC_closest(X,xi);indx2=funC_closest(X,xf);
else
    indx1=1;indx2=length(X);
end

if setting.ylim.check==1
indy1=funC_closest(Y,yi);indy2=funC_closest(Y,yf);
else
indy1=1;indy2=length(Y);
end

if setting.tlim.check==1
indt1=funC_closest(T,ti);indt2=funC_closest(T,tf);   
else
indt1=1;indt2=length(T);    
end
Xnow=X(indx1:stepS:indx2);
Ynow=Y(indy1:stepS:indy2);
Tnow=T(indt1:stepT:indt2);
[XX,YY]=meshgrid(Xnow,Ynow);

xlimm=[Xnow(1) Xnow(end)];
ylimm=[Ynow(1) Ynow(end)];
    
varb=output.eta;
if input.wall.option==1
    varb(dom.wall.charAll<1)=NaN;
end

varnow=varb(indt1:stepT:indt2,indy1:stepS:indy2,indx1:stepS:indx2);
set(h.monitorbox,'string','');
if varId==1
    var_xy=squeeze(4*sqrt(var(varnow,0,1)));
    varname='Significant wave height';
    MeanHs=mean(mean(var_xy));
    set(h.monitorbox,'string',['>> Mean Hs=',num2str(MeanHs)])
    if setting.savedata==1
       save([Proj.workdir,'PP_density_stat_data_mean_Hs.txt'],'MeanHs','-ascii');
    end
elseif varId==2
    varname='Peak period (Tp)';
    Ny=length(Ynow);Nx=length(Xnow);
    var_xy=zeros(Ny,Nx);
    
    [jProgressBar,statusbarObj]=funGui_JavProgressBar(h);
    set(jProgressBar,'Maximum',Ny*Nx, 'Value',0);
    jProgressBar.setStringPainted( true );
    statusbarObj.setText('Calculating...');
    [jbStop]=funGui_Java_stopbutton(statusbarObj);
    jProgressBar.setStringPainted( true );
    tic;
    iter=1;
    for ii=1:Ny
        for jj=1:Nx
            IDstop=funGui_eventLoopStop(jbStop);
            if IDstop==1, break;end
            
         varnow_t=squeeze(squeeze(varnow(:,ii,jj)));
         Tp=funSP_find_peak_period_Tp(Tnow,varnow_t,1);
         var_xy(ii,jj)=Tp;
         
         if mod(iter,floor(0.1*Ny*Nx))==0 || iter==floor(0.01*Ny*Nx)
             set(jProgressBar,'Maximum',Ny*Nx, 'Value',iter);
             ETA=funGui_remain_time(iter,Ny*Nx);
             statusbarObj.setText(['time remaining=', num2str(ETA)]);
         end
         iter=iter+1;
        end
    end
    set(jProgressBar,'Maximum',Ny*Nx, 'Value',Ny*Nx);
    jProgressBar.setVisible(0);
    jbStop.setVisible(0);
    statusbarObj.setText('done.');

elseif varId==3
    varname='Mean Period (Tm01)';
    Ny=length(Ynow);Nx=length(Xnow);
    var_xy=zeros(Ny,Nx);
    
    [jProgressBar,statusbarObj]=funGui_JavProgressBar(h);
    set(jProgressBar,'Maximum',Ny*Nx, 'Value',0);
    jProgressBar.setStringPainted( true );
    statusbarObj.setText('Calculating...');
    [jbStop]=funGui_Java_stopbutton(statusbarObj);
    jProgressBar.setStringPainted( true );
    tic;
    iter=1;
    for ii=1:Ny
        for jj=1:Nx
            IDstop=funGui_eventLoopStop(jbStop);
            if IDstop==1, break;end
            
         varnow_t=squeeze(squeeze(varnow(:,ii,jj)));
         Tm01=funSP_find_mean_period_Tm01(Tnow,varnow_t);
         var_xy(ii,jj)=Tm01;
         
         if mod(iter,floor(0.1*Ny*Nx))==0 || iter==floor(0.01*Ny*Nx)
             set(jProgressBar,'Maximum',Ny*Nx, 'Value',iter);
             ETA=funGui_remain_time(iter,Ny*Nx);
             statusbarObj.setText(['time remaining=', num2str(ETA)]);
         end
         iter=iter+1;
        end
    end
    set(jProgressBar,'Maximum',Ny*Nx, 'Value',Ny*Nx);
    jProgressBar.setVisible(0);
    jbStop.setVisible(0);
    statusbarObj.setText('done.');              
    
elseif varId==4
    var_xy=squeeze(mean(varnow,1));
    varname='Mean water level';    
elseif varId==5
    var_xy=squeeze(max(varnow,[],1));
    varname='Max. temporal crest';
    Etaa=output.eta(indt1:stepT:indt2,indy1:indy2,indx1:indx2);
    [valmax,Indmax]=max(reshape(Etaa,[],1));
    [indtmax,indymax,indxmax] = ind2sub(size(Etaa),Indmax);
    xx=X(indx1:indx2);yy=Y(indy1:indy2);tt=output.time(indt1:stepT:indt2);
    set(h.monitorbox,'string',['>> max elevation=',num2str(valmax),'[m], at (x,y,t)= (',num2str(xx(indxmax)),';',...
        num2str(yy(indymax)),';',num2str(tt(indtmax)),')'])
    if setting.savedata==1
        datt=[valmax xx(indxmax) yy(indymax) tt(indtmax)];
        save([Proj.workdir,'PP_density_stat_data_maxeta_xyt.txt'],'datt','-ascii');
    end
elseif varId==6
    var_xy=squeeze(min(varnow,[],1));
    varname='Min. temporal trough';
    Etaa=output.eta(indt1:stepT:indt2,indy1:indy2,indx1:indx2);
    [valmin,Indmin]=min(reshape(Etaa,[],1));
    [indtmin,indymin,indxmin] = ind2sub(size(Etaa),Indmin);
    xx=X(indx1:indx2);yy=Y(indy1:indy2);tt=output.time(indt1:stepT:indt2);
    set(h.monitorbox,'string',['>> min elevation=',num2str(valmin),'[m], at (x,y,t)= (',num2str(xx(indxmin)),';',...
        num2str(yy(indymin)),';',num2str(tt(indtmin)),')'])
    if setting.savedata==1
        datt=[valmin xx(indxmin) yy(indymin) tt(indtmin)];
        save([Proj.workdir,'PP_density_stat_data_mineta_xyt.txt'],'datt','-ascii');
    end
elseif varId==7
    varname='Average Temporal Crest (ATC)';
    Ny=length(Ynow);Nx=length(Xnow);
    var_xy=zeros(Ny,Nx);
    
    [jProgressBar,statusbarObj]=funGui_JavProgressBar(h);
    set(jProgressBar,'Maximum',Ny*Nx, 'Value',0);
    jProgressBar.setStringPainted( true );
    statusbarObj.setText('Calculating...');
    [jbStop]=funGui_Java_stopbutton(statusbarObj);
    jProgressBar.setStringPainted( true );
    tic;
    iter=1;
    for ii=1:Ny
        for jj=1:Nx
            IDstop=funGui_eventLoopStop(jbStop);
            if IDstop==1, break;end
            
         varnow_t=squeeze(squeeze(varnow(:,ii,jj)));
         [ATC,~]=funSP_find_average_temporal_crest_and_trough(Tnow,varnow_t);
         var_xy(ii,jj)=ATC;
         
         if mod(iter,floor(0.1*Ny*Nx))==0 || iter==floor(0.01*Ny*Nx)
             set(jProgressBar,'Maximum',Ny*Nx, 'Value',iter);
             ETA=funGui_remain_time(iter,Ny*Nx);
             statusbarObj.setText(['time remaining=', num2str(ETA)]);
         end
         iter=iter+1;
        end
    end
    set(jProgressBar,'Maximum',Ny*Nx, 'Value',Ny*Nx);
    statusbarObj.setText('');
    jProgressBar.setVisible(0);
    jbStop.setVisible(0);
elseif varId==8
    varname='Average Temporal Trough (ATT)';
    Ny=length(Ynow);Nx=length(Xnow);
    var_xy=zeros(Ny,Nx);
    
    [jProgressBar,statusbarObj]=funGui_JavProgressBar(h);
    set(jProgressBar,'Maximum',Ny*Nx, 'Value',0);
    jProgressBar.setStringPainted( true );
    statusbarObj.setText('Calculating...');
    [jbStop]=funGui_Java_stopbutton(statusbarObj);
    jProgressBar.setStringPainted( true );
    tic;
    iter=1;
    for ii=1:Ny
        for jj=1:Nx
            IDstop=funGui_eventLoopStop(jbStop);
            if IDstop==1, break;end
            
         varnow_t=squeeze(squeeze(varnow(:,ii,jj)));
         [~,ATT]=funSP_find_average_temporal_crest_and_trough(Tnow,varnow_t);
         var_xy(ii,jj)=ATT;
         
         if mod(iter,floor(0.1*Ny*Nx))==0 || iter==floor(0.01*Ny*Nx)
             set(jProgressBar,'Maximum',Ny*Nx, 'Value',iter);
             ETA=funGui_remain_time(iter,Ny*Nx);
             statusbarObj.setText(['time remaining=', num2str(ETA)]);
         end
         iter=iter+1;
        end
    end
    set(jProgressBar,'Maximum',Ny*Nx, 'Value',Ny*Nx);
    statusbarObj.setText('');
    jProgressBar.setVisible(0);
    jbStop.setVisible(0);
end


varlabel=[varname,' [m]'];
% if isfield(simuldata,'bdyassim')
% if input.bdyassim.option==1
% %var_xy(bdyassim.charupdate>0)=NaN;
% end
% end

if strcmpi(input.wall.option,'Yes')
    var_xy(dom.wall.charAll(indy1:stepS:indy2,indx1:stepS:indx2)<1)=NaN;
end
if setting.clim.check==0
    maxvar=max(max(max(var_xy)));
    minvar=min(min(min(var_xy)));
    climm=[minvar maxvar];
else
    climm=setting.clim.param';
end
cla(axesfig);
set(hfig,'Renderer','zbuffer'); %due to graphics driver
surf(axesfig,XX,YY,var_xy,'edgecolor','none');
xlabel(axesfig,'x[m]');ylabel(axesfig,'y[m]');zlabel(axesfig,varlabel);
xlim(axesfig,xlimm);ylim(axesfig,ylimm);
grid(axesfig,'off')
try
    caxis(axesfig,climm);
catch
    
end

try
    view(axesfig,viewplot);
catch
    view(axesfig,2);
end
colormap(axesfig,cell2mat(setting.colormap));
axes_properties(axesfig,1);
axes(axesfig)
cb=colorbar;
ylabel(cb,varlabel,'fontweight', 'bold')
axes_properties_cb(cb,1);

if setting.savedata==1
    xx=Xnow;yy=Ynow;
    varr=zeros(length(yy)+1,length(xx)+1);
    varr(1,2:end)=xx;varr(2:end,1)=yy;
    varr(2:end,2:end)=var_xy;
    save([Proj.workdir,'PP_density_stat_data_',varname,'.txt'],'varr','-ascii');
end

if CheckId.SaveFig==1
    if strcmp(CheckId.SaveFig_type,'.eps')
        ax_old = axesfig;
        f_new = figure('unit','normalized','position',[1.1 1.1 0.7 0.7]);
        set(f_new ,'Renderer','zbuffer'); %due to graphics driver
        ax_new = copyobj(ax_old,f_new);
        set(f_new,'visible','off')
        colormap(ax_new,cell2mat(setting.colormap));
        axes(ax_new)
        cb=colorbar;
        ylabel(cb,varname,'fontweight', 'bold')
        axes_properties(ax_new,1);
        axes_properties_cb(cb,1);
        set(f_new,'unit','normalized','position',[1.1 1.1 0.7 0.7])
        set(ax_new,'outerposition',[0 0 1 1]);
        saveas(ax_new,cell2mat([Proj.workdir,'PP_density_stat_',varname,CheckId.SaveFig_type]),'epsc');
        close
    else
        ax_old = axesfig;
        f_new = figure('unit','normalized','position',[1.1 1.1 0.7 0.7]);
        set(f_new ,'Renderer','zbuffer'); %due to graphics driver
        ax_new = copyobj(ax_old,f_new);
        set(f_new,'visible','off')
        colormap(ax_new,cell2mat(setting.colormap));
        axes(ax_new)
        cb=colorbar;
        ylabel(cb,varname,'fontweight', 'bold')
        axes_properties(ax_new,1);
        axes_properties_cb(cb,1);
        set(f_new,'unit','normalized','position',[1.1 1.1 0.7 0.7])
        set(ax_new,'outerposition',[0 0 1 1]);
        saveas(ax_new,cell2mat([Proj.workdir,'PP_density_stat_',varname,CheckId.SaveFig_type]))
        close;
    end
end
end


