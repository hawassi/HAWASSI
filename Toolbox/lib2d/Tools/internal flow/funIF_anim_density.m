function funIF_anim_density(h,inputvar,propdata,setting)
hfig=h.fig;
axesfig=h.IF_anim_density_axes;
set(h.monitorbox,'foregroundcolor','k','string',['>> processing...']);

if setting.view.check==1
    viewplot=setting.view.param';
else,    viewplot=2;    end

stepS=1;

if setting.coarse.check==0
    stepT=1;
else, stepT=setting.coarse.param; end

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
if setting.tlim.check==1
    indt1=funC_closest(propdata.time,setting.tlim.param(1));
    indt2=funC_closest(propdata.time,setting.tlim.param(2));
    Tnow=propdata.time(indt1:stepT:indt2);
else
    Tnow=propdata.time;
end

if setting.saveanim.check==1
    GIF_delaytime=setting.gifset.param(1);
    GIF_loopcount=setting.gifset.param(2);
    saveanim=1;
else
    saveanim=0;
end

if inputvar.axesoptId==1
        indz=funC_closest(propdata.z,inputvar.var2);
        xlab='x[m]'; ylab= 'y[m]';
        [XX,YY]=meshgrid(propdata.x(1:stepS:end),propdata.y(1:stepS:end));
    elseif inputvar.axesoptId==2
        indy=funC_closest(propdata.y,inputvar.var2);
        xlab='x[m]'; ylab= 'z[m]';
        [XX,YY]=meshgrid(propdata.x(1:stepS:end),propdata.z(1:stepS:end));
    elseif inputvar.axesoptId==3
        indx=funC_closest(propdata.x,inputvar.var2);
         xlab='y[m]'; ylab= 'z[m]';
        [XX,YY]=meshgrid(propdata.y(1:stepS:end),propdata.z(1:stepS:end));
end
    
cla(axesfig);
axes(axesfig)
if saveanim==1
    ax = axesfig;
    set(ax,'units','pixels');
    pos =get(ax,'Position');
    ti = get(ax,'TightInset');
    marg=ti(1);
    rect=[-marg, -marg, pos(3)+4*marg, pos(4)+1.5*marg];
    set(ax,'units','normalized');
end

 varlabel=StringVar{IdVar};
for ii=1:length(Tnow)
    indt=funC_closest(propdata.time,Tnow(ii));
    
    if inputvar.axesoptId==1
        PropValSqu=squeeze(squeeze(PropVal(indt,indz,1:stepS:end,1:stepS:end)));
    elseif inputvar.axesoptId==2
        PropValSqu=squeeze(squeeze(PropVal(indt,1:stepS:end,indy,1:stepS:end)));
    elseif inputvar.axesoptId==3
        PropValSqu=squeeze(squeeze(PropVal(indt,1:stepS:end,1:stepS:end,indx)));
    end
    
    
    set(hfig,'Renderer','zbuffer'); %due to graphics driver
    surf(axesfig,XX,YY,PropValSqu,'edgecolor','none');
    stringtitle=['Evolution of ',varlabel,' @ time: ',num2str(roundn(Tnow(ii),-2)), ' [s]'];
    title(axesfig,stringtitle);
    
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

      
    colormap(axesfig,cell2mat(setting.colormap));
    axes_properties(axesfig,1);
    cb=colorbar;
    ylabel(cb,varlabel,'fontweight', 'bold')
    axes_properties_cb(cb,1);
    
    
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
    
    drawnow limitrate;
    if saveanim==1
        drawnow limitrate;
        frame = getframe(ax,rect);
        pause(GIF_delaytime)
        
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        
        filename=[setting.workdir,'IF_anim_density_',varname,'.gif'];
        if ii == 1
            imwrite(imind,cm,filename,'gif','DelayTime',GIF_delaytime,'loopcount',GIF_loopcount);
        else
            imwrite(imind,cm,filename,'gif','DelayTime',GIF_delaytime, 'Writemode', 'append');
        end
    else
        drawnow limitrate;
    end
    
    
    
    if get(h.IF_anim_pause,'userdata')==1
        uiwait(gcf);
    end
    
    if get(h.IF_anim_stop,'userdata')==1
        set(h.IF_anim_stop,'userdata',0);
        break;
    end
end
set(h.monitorbox,'foregroundcolor','k','string',['>> ']);
