function fun_plotting_density_IF2D(h,inputvar,propdata,setting)
hfig=h.fig;
axesfig=h.IF_plot_density_profile_axes;
set(h.monitorbox,'foregroundcolor','k','string',['>> processing...']);



if setting.view.check==1
viewplot=setting.view.param';
else,    viewplot=2;    end

if setting.coarse.check==0
stepS=1;
else, stepS=setting.coarse.param; end


StringVar={'Velocity in x direction [m/s]','Velocity in y direction [m/s]','Velocity in z direction [m/s]','|Velocity| [m/s]' ...
    'Acceleration in x direction [m/s^2]', 'Acceleration in y direction [m/s^2]', 'Acceleration in z direction [m/s^2]', 'Acceleration [m/s^2]',...
    'Pressure (dyn.) [N/m^2]', 'Pressure (Tot.) [N/m^2]' };

if inputvar.propertiesId==1
    PropVal=propdata.dxPhi;IdVar=1;
    varname='HorzVelX';
elseif inputvar.propertiesId==2
    PropVal=propdata.dyPhi;IdVar=2;
    varname='HorzVelY';
elseif inputvar.propertiesId==3
    PropVal=propdata.dzPhi;IdVar=3;
    varname='VertVel';
elseif inputvar.propertiesId==4
    PropVal =sqrt(propdata.dxPhi.^2+propdata.dyPhi.^2+propdata.dzPhi.^2);IdVar=4;
    varname='absVel';
elseif inputvar.propertiesId==5
    PropVal=propdata.dt_dxPhi;IdVar=5;
    varname='HorzAccelX';
elseif inputvar.propertiesId==6
    PropVal=propdata.dt_dyPhi;IdVar=6;
     varname='HorzAccelY';
elseif inputvar.propertiesId==7
    PropVal=propdata.dt_dzPhi;IdVar=7;
    varname='VertAccel';
elseif inputvar.propertiesId==8
    PropVal=sqrt(propdata.dt_dxPhi.^2+propdata.dt_dyPhi.^2+propdata.dt_dzPhi.^2);IdVar=8;
    varname='absAccel';
elseif inputvar.propertiesId==9
    PropVal=propdata.dtPhi;IdVar=9;
    varname='dynPressure';
elseif inputvar.propertiesId==10
    PropVal=propdata.Ptot;IdVar=10;
     varname='totPressure';
end

if inputvar.axesoptId==1
 indt=funC_closest(propdata.time,inputvar.var1);
 indz=funC_closest(propdata.z,inputvar.var2);
 PropValSqu=squeeze(squeeze(PropVal(indt,indz,1:stepS:end,1:stepS:end)));
 xlab='x[m]'; ylab= 'y[m]';
 [XX,YY]=meshgrid(propdata.x(1:stepS:end),propdata.y(1:stepS:end));
 snapVartext1='_t'; snapVartext2='z';
elseif inputvar.axesoptId==2
 indt=funC_closest(propdata.time,inputvar.var1);
 indy=funC_closest(propdata.y,inputvar.var2);
 PropValSqu=squeeze(squeeze(PropVal(indt,1:stepS:end,indy,1:stepS:end)));  
 xlab='x[m]'; ylab= 'z[m]';
 [XX,YY]=meshgrid(propdata.x(1:stepS:end),propdata.z(1:stepS:end));
 snapVartext1='_t'; snapVartext2='y';
elseif inputvar.axesoptId==3
 indt=funC_closest(propdata.time,inputvar.var1);
 indx=funC_closest(propdata.x,inputvar.var2);
 PropValSqu=squeeze(squeeze(PropVal(indt,1:stepS:end,1:stepS:end,indx)));
 xlab='y[m]'; ylab= 'z[m]';
 [XX,YY]=meshgrid(propdata.y(1:stepS:end),propdata.z(1:stepS:end));
 snapVartext1='_t'; snapVartext2='x';
 elseif inputvar.axesoptId==4
 indy=funC_closest(propdata.y,inputvar.var1);
 indz=funC_closest(propdata.z,inputvar.var2);
 PropValSqu=squeeze(squeeze(PropVal(1:stepS:end,indz,indy,1:stepS:end))); 
 xlab='x[m]'; ylab= 't[m]';
 [XX,YY]=meshgrid(propdata.x(1:stepS:end),propdata.time(1:stepS:end));
 snapVartext1='_y'; snapVartext2='z';
 elseif inputvar.axesoptId==5
 indx=funC_closest(propdata.x,inputvar.var1);
 indz=funC_closest(propdata.z,inputvar.var2);
 PropValSqu=squeeze(squeeze(PropVal(1:stepS:end,indz,1:stepS:end,indx))).'; 
 xlab='t[m]'; ylab= 'y[m]';
 [XX,YY]=meshgrid(propdata.time(1:stepS:end),propdata.y(1:stepS:end));
 snapVartext1='_x'; snapVartext2='z';
 elseif inputvar.axesoptId==6
 indx=funC_closest(propdata.x,inputvar.var1);
 indy=funC_closest(propdata.y,inputvar.var2);
 PropValSqu=squeeze(squeeze(PropVal(1:stepS:end,1:stepS:end,indy,indx))).'; 
 xlab='t[m]'; ylab= 'z[m]';
 [XX,YY]=meshgrid(propdata.time(1:stepS:end),propdata.z(1:stepS:end));
 snapVartext1='_x'; snapVartext2='y';
end

cla(axesfig);
set(hfig,'Renderer','zbuffer'); %due to graphics driver
surf(axesfig,XX,YY,PropValSqu,'edgecolor','none');
xlabel(xlab);ylabel(ylab);
grid(axesfig,'off')

try
view(axesfig,viewplot);
catch
view(axesfig,2);    
end


if setting.clim.check==1
    climm=setting.clim.param';
    zlim(axesfig,climm);caxis(axesfig,climm);
end

if setting.xlim.check==1
xi=setting.xlim.param(1);xf=setting.xlim.param(2);
xlim(axesfig,[xi xf]);
end
if setting.ylim.check==1
yi=setting.ylim.param(1);yf=setting.ylim.param(2);
ylim(axesfig,[yi yf]);
end
varlabel=StringVar{IdVar};
colormap(axesfig,cell2mat(setting.colormap));
axes_properties(axesfig,1);
axes(axesfig)
cb=colorbar;
ylabel(cb,varlabel,'fontweight', 'bold')
axes_properties_cb(cb,1);
set(h.monitorbox,'foregroundcolor','k','string',['>> ']);


if inputvar.levelCb==1
    if setting.level.check==1
        thresholdLevel=setting.level.param;
    else
        thresholdLevel=[];
    end
    hold(axesfig,'on');
%     if length(thresholdLevel)>10
%         thresholdLevel= thresholdLevel(1:2:end);
%     end
    funP_contour(axesfig,XX,YY,PropValSqu,thresholdLevel,0,'k',1.5);
    hold(axesfig,'off');
end


if setting.savefig.check==1
    if strcmp(setting.savefig.format,'.eps')
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
        saveas(ax_new,cell2mat([setting.workdir,'IF_density_profile_',varname,snapVartext1,num2str(inputvar.var1),'_',snapVartext2,num2str(inputvar.var2),CheckId.SaveFig_type]),'epsc');
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
        saveas(ax_new,cell2mat([setting.workdir,'IF_density_profile_',varname,snapVartext1,num2str(inputvar.var1),'_',snapVartext2,num2str(inputvar.var2),setting.savefig.format]))
        close(f_new);
    end
end
set(h.monitorbox,'foregroundcolor','k','string',['>>']);