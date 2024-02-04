function funP_line_plot_statistic(spatsnap,simuldata,setting,axesfig,h)
varId=setting.var;
output=simuldata.output;
input=simuldata.input;
dom=simuldata.dom;
Proj=simuldata.Proj;
X=output.X;Y=output.Y;
T=output.time;
if setting.coarse.check==0
stepS=1;
else
stepS=setting.coarse.param;  
end
CheckId.SaveFig=setting.savefig.check;
CheckId.SaveFig_type=setting.savefig.format;



if spatsnap.y_check==1
if setting.spatlim.check==1   
spati=setting.spatlim.param(1);
spatf=setting.spatlim.param(2);    
indspat1=funC_closest(X,spati);indspat2=funC_closest(X,spatf);
else
indspat1=1; indspat2=length(X);   
end
spat1=X;spat2=Y;spat2snap=spatsnap.y_snap;
else
if setting.spatlim.check==1   
spati=setting.spatlim.param(1);
spatf=setting.spatlim.param(2);    
indspat1=funC_closest(Y,spati);indspat2=funC_closest(Y,spatf);
else
indspat1=1; indspat2=length(Y);       
end
spat1=Y;spat2=X;spat2snap=spatsnap.x_snap;
end

if setting.spatlim.check==1
xlimm=[spat1(indspat1) spat1(indspat2)];
else
xlimm=[spat1(1) spat1(end)];   
end
spat=spat1(indspat1:stepS:indspat2);

if setting.tlim.check==1
indt1=funC_closest(T,setting.tlim.param(1));
indt2=funC_closest(T,setting.tlim.param(2));
else
indt1=1;indt2=length(T);   
end

cla(axesfig);
var_xy=zeros(size(spat));
Nspat=length(spat);
[jProgressBar,statusbarObj]=funGui_JavProgressBar(h);
set(jProgressBar,'Maximum',Nspat, 'Value',0);
jProgressBar.setStringPainted( true );
statusbarObj.setText('Plotting...');
[jbStop]=funGui_Java_stopbutton(statusbarObj);
jProgressBar.setStringPainted( true );
tic;
iter=1;

for ii=1:Nspat
    IDstop=funGui_eventLoopStop(jbStop);
    if IDstop==1, break;end
    
    if spatsnap.y_check==1
        indx=funC_closest(X,spat(ii));
        indy=funC_closest(Y,spat2snap);
    else
        indx=funC_closest(X,spat2snap);
        indy=funC_closest(Y,spat(ii));
    end
    
    if varId==1
        var_xy(ii)=squeeze(4*sqrt(var(output.eta(indt1:indt2,indy,indx),0,1)));
        varname='Significant wave height';
        varlabel=[varname, ' [m]'];
    elseif varId==2
        var_xy(ii)=squeeze(max(output.eta(indt1:indt2,indy,indx),[],1));
        varname='Maximum temporal crest';  varlabel=[varname, ' [m]'];
    elseif varId==3
        var_xy(ii)=squeeze(min(output.eta(indt1:indt2,indy,indx),[],1));
        varname='Minimum temporal trough';  varlabel=[varname, ' [m]'];
    elseif varId==4
        varbnow_t=squeeze(squeeze(output.eta(indt1:indt2,indy,indx)));Tnow=T(indt1:indt2);
        [var_xy(ii),~]=funSP_find_average_temporal_crest_and_trough(Tnow,varbnow_t);
        varname='Average temporal crest';  varlabel=[varname, ' [m]'];
    elseif varId==5
        varbnow_t=squeeze(squeeze(output.eta(indt1:indt2,indy,indx)));Tnow=T(indt1:indt2);
        [~,var_xy(ii)]=funSP_find_average_temporal_crest_and_trough(Tnow,varbnow_t);
        varname='Average temporal trough';  varlabel=[varname, ' [m]'];
    elseif varId==6
        var_xy(ii)=squeeze(mean(output.eta(indt1:indt2,indy,indx),1));
        varname='Mean water level'; varlabel=[varname, ' [m]'];
    elseif varId==7
        etanow=squeeze(output.eta(indt1:indt2,indy,indx));
        Tnow=T(indt1:indt2);
        var_xy(ii)=funSP_find_peak_period(Tnow,etanow);
        varname='Peak period'; varlabel=[varname, ' [s]'];
    elseif varId==8
        etanow=squeeze(output.eta(indt1:indt2,indy,indx));
        var_xy(ii)=funSP_find_skewness(etanow);
        varname='Skewness'; varlabel=[varname];
    elseif varId==9
        etanow=squeeze(output.eta(indt1:indt2,indy,indx));
        var_xy(ii)=-funSP_find_asymmetry(etanow);
        varname='-Asymmetry';varlabel=[varname];
    elseif varId==10
        etanow=squeeze(output.eta(indt1:indt2,indy,indx));
        var_xy(ii)=kurtosis(etanow-mean(etanow));
        varname='Kurtosis'; varlabel=[varname];
    end
    if mod(iter,floor(0.1*Nspat))==0 || iter==floor(0.01*Nspat)
        set(jProgressBar,'Maximum',Nspat, 'Value',iter);
        ETA=funGui_remain_time(iter,Nspat);
        statusbarObj.setText(['time remaining=', num2str(ETA)]);
    end
    iter=iter+1;
end
set(jProgressBar,'Maximum',Nspat, 'Value',Nspat);
jProgressBar.setVisible(0);
jbStop.setVisible(0);
statusbarObj.setText('done.');

if spatsnap.y_check==1
xlab='x[m]';labb=['_at_y=',num2str(spat2snap),'_'];
else
xlab='y[m]';labb=['_at_x=',num2str(spat2snap),'_'];
end

plot(axesfig,spat,var_xy,'or');

if setting.savedat==1
   varr=[spat' var_xy'];
  save([Proj.workdir,'PP_line_statistic',varname,labb,'.txt'],'varr','-ascii');       
end

if setting.zlim.check==0
maxvar=max(max(max(var_xy))); 
minvar=min(min(min(var_xy))); 
zlimm=[minvar maxvar];
if minvar==maxvar
zlimm=[minvar-0.01 maxvar+0.01];    
end
else
zlimm=setting.zlim.param';    
end
axes_properties(axesfig,1);
xlabel(axesfig,xlab);ylabel(axesfig,varlabel);
xlim(axesfig,xlimm);ylim(axesfig,zlimm);
axes_properties(axesfig,1);

if CheckId.SaveFig==1
    if strcmp(CheckId.SaveFig_type,'.eps')
        ax_old = axesfig;
        f_new = figure('unit','normalized','position',[1.1 1.1 0.7 0.7]);
        ax_new = copyobj(ax_old,f_new);
        axes_properties(ax_new,1);
        set(f_new,'unit','normalized','position',[1.1 1.1 0.7 0.7])
        set(ax_new,'outerposition',[0 0 1 1]);
        saveas(ax_new,cell2mat([Proj.workdir,'PP_line_statistic',varname,labb,CheckId.SaveFig_type]),'epsc');
        close
    else
        ax_old = axesfig;
        f_new = figure('unit','normalized','position',[1.1 1.1 0.7 0.7]);
        ax_new = copyobj(ax_old,f_new);
        axes_properties(ax_new,1);
        set(f_new,'unit','normalized','position',[1.1 1.1 0.7 0.7])
        set(ax_new,'outerposition',[0 0 1 1]);
        saveas(ax_new,cell2mat([Proj.workdir,'PP_line_statistic',varname,labb,CheckId.SaveFig_type]))
        close;
    end
end

end