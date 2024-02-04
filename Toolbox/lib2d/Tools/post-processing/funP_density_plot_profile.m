function funP_density_plot_profile(tsnap,simuldata,setting,axesfig,handles)
hfig=handles.fig;

varId=setting.var;
output=simuldata.output;
input=simuldata.input;
influx=simuldata.influx;
model=simuldata.model;
bathy=simuldata.dom.bathy.profile;

if isfield(simuldata,'ivp')
ivp=simuldata.ivp;
end

if isfield(simuldata,'bdyassim')
bdyassim=simuldata.bdyassim;
end

dom=simuldata.dom;
Proj=simuldata.Proj;

X=output.X;Y=output.Y;
T=output.time;
dT=T(2)-T(1);

if setting.xlim.check==1
xi=setting.xlim.param(1);
xf=setting.xlim.param(2);
else, xi=X(1); xf=X(end); 
end

if setting.ylim.check==1
yi=setting.ylim.param(1);
yf=setting.ylim.param(2);
else, yi=Y(1); yf=Y(end); 
end


if setting.view.check==1
viewplot=setting.view.param';
else    
viewplot=2;    
end

if setting.coarse.check==0
stepS=1;
else
stepS=setting.coarse.param;  
end

CheckId.SaveFig=setting.savefig.check;
CheckId.SaveFig_type=setting.savefig.format;

indx1=funC_closest(X,xi);indx2=funC_closest(X,xf);
indy1=funC_closest(Y,yi);indy2=funC_closest(Y,yf);

[XX,YY]=meshgrid(X(indx1:stepS:indx2),Y(indy1:stepS:indy2));

xlimm=[X(indx1) X(indx2)];
ylimm=[Y(indy1) Y(indy2)];

indt=funC_closest(T,tsnap);

kx           = funC_freqspace(X)';
ky           = funC_freqspace(Y)';
[Kx,Ky]  = meshgrid(kx,ky);

if varId==1
var_xy=squeeze(output.eta(indt,:,:));
varname='elevation';
varlabel=[varname,' [m]'];
elseif varId==2
eta=squeeze(output.eta(indt,:,:));
gradeta=funOprt_grad2d(Kx,Ky,fft2(eta)); 
%var_xy=sqrt(gradeta.x.^2+gradeta.y.^2);
var_xy=gradeta.x+gradeta.y;
varname='gradient of elevation';   
varlabel=[varname];
elseif varId==3
    if model.phiForm==1
        var_xy=squeeze(output.phi(indt,:,:));
    else
        Velx=squeeze(output.u(indt,:,:));%vel.x(indy1:stepS:indy2,indx1:stepS:indx2);
        Vely=squeeze(output.v(indt,:,:));
        var_xy=intgrad2(Velx,Vely,X(2)-X(1),Y(2)-Y(1));
    end
varname='potential';  
varlabel=[varname,' [m^2/s]'];
elseif varId==4
    if model.phiForm==1
        phi=squeeze(output.phi(indt,:,:));
        vel=funOprt_grad2d(Kx,Ky,fft2(phi));
        var_xy=sqrt(vel.x.^2+vel.y.^2);
        Velx=vel.x(indy1:stepS:indy2,indx1:stepS:indx2);
        Vely=vel.y(indy1:stepS:indy2,indx1:stepS:indx2);
    else
        Velx=squeeze(output.u(indt,:,:));%vel.x(indy1:stepS:indy2,indx1:stepS:indx2);
        Vely=squeeze(output.v(indt,:,:));
        var_xy=sqrt(Velx.^2+Vely.^2);
        Velx=Velx(indy1:stepS:indy2,indx1:stepS:indx2);
        Vely=Vely(indy1:stepS:indy2,indx1:stepS:indx2);
    end
varname='|velocity|';  
varlabel=[varname,' [m/s]'];
elseif varId==5
    if model.phiForm==1
        phi=squeeze(output.phi(indt,:,:));
        vel=funOprt_grad2d(Kx,Ky,fft2(phi));
        var_xy=vel.x;
        Velx=vel.x(indy1:stepS:indy2,indx1:stepS:indx2);
        Vely=vel.y(indy1:stepS:indy2,indx1:stepS:indx2);
    else
        Velx=squeeze(output.u(indt,:,:));%vel.x(indy1:stepS:indy2,indx1:stepS:indx2);
        Vely=squeeze(output.v(indt,:,:));
        var_xy=Velx;
        Velx=Velx(indy1:stepS:indy2,indx1:stepS:indx2);
        Vely=Vely(indy1:stepS:indy2,indx1:stepS:indx2);
    end
varname='velocity (in x dir)';  
varlabel=[varname,' [m/s]'];
elseif varId==6
    if model.phiForm==1
        phi=squeeze(output.phi(indt,:,:));
        vel=funOprt_grad2d(Kx,Ky,fft2(phi));
        var_xy=vel.y;
        Velx=vel.x(indy1:stepS:indy2,indx1:stepS:indx2);
        Vely=vel.y(indy1:stepS:indy2,indx1:stepS:indx2);
    else
        Velx=squeeze(output.u(indt,:,:));%vel.x(indy1:stepS:indy2,indx1:stepS:indx2);
        Vely=squeeze(output.v(indt,:,:));
        var_xy=Vely;
        Velx=Velx(indy1:stepS:indy2,indx1:stepS:indx2);
        Vely=Vely(indy1:stepS:indy2,indx1:stepS:indx2);
    end
varname='velocity (in y dir)';  
varlabel=[varname,' [m/s]'];
end
    
% if isfield(simuldata,'bdyassim')
% if input.bdyassim.option==1
% var_xy(bdyassim.charupdate>0)=NaN;
% end
% end

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

if strcmpi(setting.colormap,'jet')
  colE='c';
else
  colE='r';  
end

if setting.level.eta_check==1
   varlevel=squeeze(output.eta(indt,:,:));
   flag_level=1;
elseif setting.level.phi_check==1
    if model.phiForm==1
        varlevel=squeeze(output.phi(indt,:,:));
    else
        Velx=squeeze(output.u(indt,:,:));%vel.x(indy1:stepS:indy2,indx1:stepS:indx2);
        Vely=squeeze(output.v(indt,:,:));
        varlevel=intgrad2(Velx,Vely,X(2)-X(1),Y(2)-Y(1));
    end
   flag_level=1;
else
   flag_level=0; 
end

if setting.level.extremecrest==1 && setting.level.extremetrough==0
    varExtreme=squeeze(output.eta(indt,:,:));
    flag_extreme=1;
elseif setting.level.extremetrough==1 && setting.level.extremecrest==0
    varExtreme=squeeze(output.eta(indt,:,:));
    flag_extreme=1;
elseif setting.level.extremecrest==1 && setting.level.extremetrough==1
    varExtreme=squeeze(output.eta(indt,:,:));
    flag_extreme=2;
else
    flag_extreme=0;
end

if strcmpi(input.wall.option,'Yes')
    var_xy(dom.wall.charAll<1)=NaN;
    if flag_level==1
    varlevel(dom.wall.charAll<1)=NaN;
    end
    if flag_extreme==1
    varExtreme(dom.wall.charAll<1)=NaN;
    end
end
var_xy=var_xy(indy1:stepS:indy2,indx1:stepS:indx2);

if flag_level==1
varlevel=varlevel(indy1:stepS:indy2,indx1:stepS:indx2);
end
if flag_extreme>0
varExtreme=varExtreme(indy1:stepS:indy2,indx1:stepS:indx2);
end



cla(axesfig);
set(hfig,'Renderer','zbuffer'); %due to graphics driver
surf(axesfig,XX,YY,var_xy,'edgecolor','none');
grid(axesfig,'off')


if setting.level.extremecrest==1  && setting.level.extremetrough==0
etaCheck=varExtreme;
IdExtreme =find(etaCheck>Hs);
colr='k';
if isempty(IdExtreme)
    flag_extreme=0;
    set(handles.monitorbox,'string','>>There is no extreme crest','foregroundcolor','k')
else
    set(handles.monitorbox,'string','>>','foregroundcolor','k');
end
elseif setting.level.extremetrough==1 && setting.level.extremecrest==0
etaCheck=varExtreme;
IdExtreme =find(etaCheck<-Hs);
colr='w';
if isempty(IdExtreme)
    flag_extreme=0;
    set(handles.monitorbox,'string','>>There is no extreme trough','foregroundcolor','k')
else
    set(handles.monitorbox,'string','>>','foregroundcolor','k');
end
elseif setting.level.extremecrest==1 && setting.level.extremetrough==1 
etaCheck=varExtreme;
IdExtremeMax=find(etaCheck>Hs);
IdExtremeMin=find(etaCheck<-Hs);
  if ~isempty(IdExtremeMax)&& isempty(IdExtremeMin)
      flag_extreme=1;
      IdExtreme=IdExtremeMax;
      colr='k';
      set(handles.monitorbox,'string','>>There is no extreme trough','foregroundcolor','k')
  elseif isempty(IdExtremeMax)&& ~isempty(IdExtremeMin)
      flag_extreme=1;
      IdExtreme=IdExtremeMin;
      colr='w';
      set(handles.monitorbox,'string','>>There is no extreme crest','foregroundcolor','k')     
  elseif ~isempty(IdExtremeMax)&& ~isempty(IdExtremeMin)
      set(handles.monitorbox,'string','>>','foregroundcolor','k');
      flag_extreme=2;
  else
     flag_extreme=0; 
     set(handles.monitorbox,'string','>>There is no extreme crest/trough','foregroundcolor','k')
  end
end

if flag_extreme==1
hold(axesfig,'on');
axes(axesfig);
plot3(axesfig,XX(IdExtreme),YY(IdExtreme),etaCheck(IdExtreme),'o',...
    'MarkerFaceColor',colr,'MarkerEdgeColor',colE);
hold(axesfig,'off');
end

if flag_extreme==2
        hold(axesfig,'on');
        plot3(axesfig,XX(IdExtremeMax),YY(IdExtremeMax),etaCheck(IdExtremeMax),'o',...
            'MarkerFaceColor','k','MarkerEdgeColor',colE);
        hold(axesfig,'off');
        hold(axesfig,'on');
        plot3(axesfig,XX(IdExtremeMin),YY(IdExtremeMin),etaCheck(IdExtremeMin),'o',...
            'MarkerFaceColor','w','MarkerEdgeColor',colE);
        hold(axesfig,'off');
end

if flag_level==1
    if setting.levelquiver.check==1
        thresholdLevel=setting.levelquiver.param;
    else        
        thresholdLevel=[];
    end
    hold(axesfig,'on');
%     if length(thresholdLevel)>10
%         thresholdLevel= thresholdLevel(1:2:end);
%     end
    funP_contour(axesfig,XX,YY,varlevel,thresholdLevel,max(max(max(var_xy))),'k',1.5);
    hold(axesfig,'off');
end

if setting.level.quiver_check==1
hold(axesfig,'on');
axes(axesfig);
if setting.levelquiver.check==1  
scale=setting.levelquiver.param; 
else
scale=1;    
end
ZZ=var_xy;
quiver3(axesfig,XX,YY,ZZ,Velx,Vely,zeros(size(ZZ)),scale,'color','k');
hold(axesfig,'off');
end

if setting.bathy.cb==1
    hold(axesfig,'on');
    Bscale=setting.bathy.scale;
    if   viewplot==2, viewplot=3; end
    plot3(axesfig,XX,YY,bathy(indy1:stepS:indy2,indx1:stepS:indx2)*Bscale,'color',[0.2 0.2 0.2])
    hold(axesfig,'off');
end


xlabel(axesfig,'x[m]');ylabel(axesfig,'y[m]');zlabel(axesfig,varlabel);

xlim(axesfig,xlimm);ylim(axesfig,ylimm);

if setting.clim.check==1
    climm=setting.clim.param';
    zlim(axesfig,climm);caxis(axesfig,climm);
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

if ~isempty(output.break_nodes)
    if varId==1
        if strcmp(model.breaking.check,'Yes') &&T(indt)>=output.break_nodes(1,1)...
                &&T(indt)<=output.break_nodes(end,1)
            indtBreak=funC_closest(output.break_nodes(:,1),T(indt));
            
            
            
            nodesBreaktemp=output.break_nodes(indtBreak,2:end);
            nodesBreak=nodesBreaktemp(nodesBreaktemp>0);
            
            if all(abs(T(indt)-output.break_nodes(:,1))>dT)
            else
                var_xyB=squeeze(output.eta(indt,:,:));
                hold(axesfig,'on');
                plot3(axesfig,dom.XX([nodesBreak]),dom.YY([nodesBreak]),var_xyB([nodesBreak]),'o',...
                    'MarkerFaceColor','w','MarkerEdgeColor','w');
                hold(axesfig,'off');
            end
        end
    end
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
        saveas(ax_new,cell2mat([Proj.workdir,'PP_density_profile_',varname,'_@t=',num2str(tsnap),CheckId.SaveFig_type]),'epsc');
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
        saveas(ax_new,cell2mat([Proj.workdir,'PP_density_profile_',varname,'_@t=',num2str(tsnap),CheckId.SaveFig_type]))
        close(f_new);
    end
end