function funP_density_plot_signal(spatsnap,simuldata,setting,axesfig,hfig)
varId=setting.var;
output=simuldata.output;
input=simuldata.input;
dom=simuldata.dom;
Proj=simuldata.Proj;
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
stepT=1;
stepS=1;
else
stepS=1;    
stepT=setting.coarse.param;
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
indt1=funC_closest(T,ti);indt2=funC_closest(T,tf);

[XX,YY]=meshgrid(T(1:stepT:end),spat1(1:stepS:end));
if setting.tlim.check==1
xlimm=[T(indt1) T(indt2)];
else
xlimm=[T(1) T(end)];    
end


if setting.spatlim.check==1
ylimm=[spat1(indspat1) spat1(indspat2)];
else
ylimm=[spat1(1) spat1(end)];   
end


indspatsnap=funC_closest(spat2,spat2snap);

if varId==1
if spatsnap.y_check==1
var_xy=squeeze(output.eta(1:stepT:end,indspatsnap,1:stepS:end))';
else
var_xy=squeeze(output.eta(1:stepT:end,1:stepS:end,indspatsnap))';    
end
varname='elevation';
varlabel=[varname,' [m]'];
elseif varId==2
    if simuldata.model.phiForm==1
        if spatsnap.y_check==1
            var_xy=squeeze(output.phi(1:stepT:end,indspatsnap,1:stepS:end)).';
        else
            var_xy=squeeze(output.phi(1:stepT:end,1:stepS:end,indspatsnap)).';
        end
        
        varname='potential';
        varlabel=[varname,' [m^2/s]'];
    else
        if spatsnap.y_check==1
        var_xy=squeeze(output.u(1:stepT:end,indspatsnap,1:stepS:end));%vel.x(indy1:stepS:indy2,indx1:stepS:indx2);
        else
        var_xy=squeeze(output.u(1:stepT:end,1:stepS:end,indspatsnap)).';%vel.x(indy1:stepS:indy2,indx1:stepS:indx2);
        end
        
        varname='velocity (in x dir.)';
        varlabel=[varname,' [m/s]'];
    end
elseif varId==3
    if spatsnap.y_check==1
        var_xy=squeeze(output.v(1:stepT:end,indspatsnap,1:stepS:end));%vel.x(indy1:stepS:indy2,indx1:stepS:indx2);
    else
        var_xy=squeeze(output.v(1:stepT:end,1:stepS:end,indspatsnap)).';%vel.x(indy1:stepS:indy2,indx1:stepS:indx2);
    end
    
    varname='velocity (in y dir.)';
    varlabel=[varname,' [m/s]'];
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
xlabel(axesfig,'time[s]');zlabel(axesfig,varlabel);

if spatsnap.x_check==1
    ylabel(axesfig,'y[m]');
else
    ylabel(axesfig,'x[m]');
end

xlim(axesfig,xlimm);ylim(axesfig,ylimm);caxis(axesfig,climm);
axes_properties(axesfig,1);

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
axes_properties(axesfig,1);

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
        if spatsnap.x_check==1
        saveas(ax_new,cell2mat([Proj.workdir,'PP1_',varname,'_@x=',num2str(spat2snap),CheckId.SaveFig_type]),'epsc')
        else
        saveas(ax_new,cell2mat([Proj.workdir,'PP1_',varname,'_@y=',num2str(spat2snap),CheckId.SaveFig_type]),'epsc')    
        end
        close;
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
        if spatsnap.x_check==1
        saveas(ax_new,cell2mat([Proj.workdir,'PP1_density_signal_',varname,'_@x=',num2str(spat2snap),CheckId.SaveFig_type]))
        else
        saveas(ax_new,cell2mat([Proj.workdir,'PP_density_signal_',varname,'_@y=',num2str(spat2snap),CheckId.SaveFig_type]))    
        end
        close;
    end
end
    



end