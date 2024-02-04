% !--------------------------------------------------------------------------------------
% !
% !    Copyright (C) 2024 - LabMath-Indonesia
% !
% !    This program is free software: you can redistribute it and/or modify
% !    it under the terms of the GNU General Public License as published by
% !    the Free Software Foundation, either version 3 of the License, or
% !    (at your option) any later version.
% !
% !    This program is distributed in the hope that it will be useful,
% !    but WITHOUT ANY WARRANTY; without even the implied warranty of
% !    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% !    GNU General Public License for more details.
% !
% !    You should have received a copy of the GNU General Public License
% !    along with this program.  If not, see <http://www.gnu.org/licenses/>.
% !
% !   Contributors list:
% !   - R. Kurnia
% !   - E. van Groesen
% !--------------------------------------------------------------------------------------

function PreProcSpatialView(h,preproc)

dom=preproc(1).dom;
influx=preproc(1).influx;
proj=preproc(1).Proj;
bdyassim=preproc(1).bdyassim;
input=preproc(1).input;

prevSaveDir=[proj.workdir,'\preview\'];
if ~isdir(prevSaveDir), mkdir(prevSaveDir); end
%%%%%%%%%%%%%%dispersion plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Oprt_Om2d_fun= str2func(['Oprt_',model.dispersion]);
% omAdd        = model.OmFun;



%%%%%%Spatial
% delete(findall(h.panel.prev_spatial,'type','axes'));
% h.preview.spatial_axes1= axes('parent', h.panel.prev_spatial, 'Position', [0.12 0.12 0.38 0.7]);
% %grid(h.preview.spatial_axes1,'On');
% h.preview.spatial_axes2= axes('parent', h.panel.prev_spatial, 'Position', [0.57 0.12 0.33 0.7]);
% %grid(h.preview.spatial_axes2,'On');

if log2(dom.Nx)>9,stepdx=1;else stepdx=1;end
if log2(dom.Ny)>9,stepdy=1;else stepdy=1;end

if strcmpi(input.wave.option,'Yes')
    Gam=0;
    for I=1:input.wave.N
        GamI=influx.gen(I).Gam;
        Gam=Gam+GamI./max(max(GamI(1:stepdy:end,(1:stepdx:end))));
    end
else
    Gam=zeros(size(dom.XX));
end
 SpatDom=(1-dom.cfSA)+Gam;%

% XX=reshape(dom.XX,[],1);YY=reshape(dom.YY,[],1);
% Fspatdom=scatteredInterpolant(XX,YY,reshape(SpatDom,[],1));
% x0=linspace(min(XX),max(XX),100);
% y0=linspace(min(YY),max(YY),100);
% [xx,yy]=meshgrid(x0,y0);
% xres=reshape(xx,length(y0),[]);
% yres=reshape(yy,length(y0),[]);
%Fspatdom(xres,yres);


[xx,yy]=meshgrid(dom.X(1:stepdx:end),dom.Y(1:stepdy:end));
axis1=h.preview.spatial_axes1;
SpatPlotId=get(h.preview.spatial_popup_var,'value');
hold(axis1,'off');

if strcmpi(input.wall.option,'Yes')
    WallDom=1-dom.wall.charAll;
end


if SpatPlotId==1
    surf(axis1,xx,yy,SpatDom(1:stepdy:end,1:stepdx:end),'edgecolor','none')
    title(axis1,'Fourier Boundary & Influx char.')
    if input.bdyassim.option==1
        hold(axis1,'on');
        charupdate=bdyassim.charupdate.*dom.cfSA;
        surf(axis1,xx,yy,charupdate(1:stepdy:end,1:stepdx:end),'edgecolor','none')
        hold(axis1,'off');
    end
    SavefigPlotName='Fourier_bdy_influx_char.png';
elseif SpatPlotId==2
    if max(max(Gam))~=0
    Gam=Gam./max(max(Gam));
    else
    Gam=0;
    end
    Spatplot=Gam+influx.ChiAdj;
    surf(axis1,xx,yy,Spatplot(1:stepdy:end,1:stepdx:end),'edgecolor','none')
    title(axis1,'Influx char. & nonlinear adjustment')
    if input.bdyassim.option==1
        hold(axis1,'on');
        charupdate=bdyassim.charupdate.*dom.cfSA;
        surf(axis1,xx,yy,charupdate(1:stepdy:end,1:stepdx:end),'edgecolor','none')
        hold(axis1,'off');
    end
    SavefigPlotName='influx_char_nonlinadj.png';
    grid(axis1,'Off');
elseif SpatPlotId==3
    Bathyy=dom.bathy.profile;
    if strcmpi(input.wall.option,'Yes')
        Bathyy(WallDom>0)=NaN;
    end
    surf(axis1,xx,yy,Bathyy(1:stepdy:end,1:stepdx:end),'edgecolor','none')
    zlabel(axis1,'bathymetry [m]')
    title(axis1,'Bathymetry');
    cbname='bathymetry [m]';
    SavefigPlotName='bathymetry.png';
elseif SpatPlotId==4 && get(h.bathymetry_checkbox_friction,'value')==1

    surf(axis1,xx,yy,dom.friction.Char(1:stepdy:end,1:stepdx:end),'edgecolor','none')
    zlabel(axis1,'Friction Coef')
    title(axis1,'Bottom friction coefficient');
    cbname='Friction Coef.';
    SavefigPlotName='bottom_friction.png';
elseif SpatPlotId==5 && strcmpi(input.wall.option,'Yes')
    surf(axis1,xx,yy,SpatDom(1:stepdy:end,1:stepdx:end),'edgecolor','none')
    hold(axis1,'on');
    surf(axis1,xx,yy,WallDom(1:stepdy:end,1:stepdx:end),'edgecolor','none')
    hold(axis1,'off');
    title(axis1,'Wall, Fourier bdy and Influxing characteristic function');
    SavefigPlotName='Wall_Fourierbdy_influx_char.png';
elseif SpatPlotId==6 && strcmpi(input.wall.option,'Yes') && dom.wall.NInfl>0
    surf(axis1,xx,yy,dom.wall.SkewGam.All(1:stepdy:end,1:stepdx:end),'edgecolor','none')
    title(axis1,'Wall influxing source function');
    SavefigPlotName='Wall_infl_source_char.png';
elseif SpatPlotId==7 && strcmpi(input.wall.option,'Yes') && dom.wall.NInfl>0
    surf(axis1,xx,yy,influx.ChiAdjWall(1:stepdy:end,1:stepdx:end),'edgecolor','none')
    SavefigPlotName='Wall_infl_nonlin_adj.png';
    title(axis1,'Wall influxing nonlinear adjustment');
else
   SpatPlotId=8; 
end
grid(axis1,'Off');
if SpatPlotId<8
xlim(axis1,[min(min(xx)) max(max(xx))]);
ylim(axis1,[min(min(yy)) max(max(yy))]);
xlabel(axis1,'x [m]');ylabel(axis1,'y[m]');
view(axis1,2)
axes_properties(axis1,1.5);
colormap(axis1,'jet');%colorbar;
if SpatPlotId==3 || SpatPlotId==4
    axes(axis1)
    cb=colorbar;
    ylabel(cb,cbname)
end

ax_old = axis1;
f_new = figure('position',[1.5 1.5 0.7 0.7]);
set(f_new,'Renderer','zbuffer');
ax_new = copyobj(ax_old,f_new);
set(f_new,'visible','off')
colormap(ax_new,'jet');%colorbar;
axes_properties(ax_new,1);
if SpatPlotId==3 || SpatPlotId==4
    axes(ax_new)
    cb=colorbar;ylabel(cb,cbname)
    axes_properties_cb(cb,1);
end
set(f_new,'unit','normalized','position',[1.5 1.5 0.7 0.7])
set(ax_new,'outerposition',[0 0 1 1]);
saveas(ax_new,[prevSaveDir,SavefigPlotName]);
%saveas(ax_new,[prevSaveDir,'simul_domain_char.fig']);
close(f_new)
axes_properties(axis1,1.5);
end


end