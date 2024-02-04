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

function PreProcView(h,preproc)

dom=preproc(1).dom;model=preproc(1).model;
par=preproc(1).par;influx=preproc(1).influx;
proj=preproc(1).Proj;g=par.g;
ivp=preproc(1).ivp;bdyassim=preproc(1).bdyassim;
input=preproc(1).input;log=preproc(1).log;
meandepth=mean(mean(-dom.bathy.profile));

prevSaveDir=[proj.workdir,'\preview\'];
if ~isdir(prevSaveDir), mkdir(prevSaveDir); end
%%%%%%%%%%%%%%dispersion plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Oprt_Om2d_fun= str2func(['Oprt_',model.dispersion]);
% omAdd        = model.OmFun;

Om_fun       = str2func(['funOprt_Om',model.dispersion(5:end)]);
invOm_fun    = str2func(['funOprt_invOm',model.dispersion(5:end)]);
Up_fun       = @(k,d,g,omAdd) Om_fun(k,d,g,omAdd)./k;
Ug_fun       = str2func(['funOprt_Ug',model.dispersion(5:end)]);


%Oprt_Up2d_fun= str2func(['Oprt_',model.phasevel]);

kx=linspace(0,2*max(dom.kx)/dom.cutfracx,100);%% max(dom.kx)is kmax/2
ky=linspace(0,2*max(dom.ky)/dom.cutfracy,100);
[Kx,Ky]=meshgrid(kx,ky);
KK=sqrt(Kx.^2+Ky.^2);
kk=linspace(0,max(max(KK)),100);
%Omega2D=Oprt_Om2d_fun(Kx,Ky,meandepth,g,[1 1]);
OmegaAbsK=Om_fun(kk,meandepth,g,'');
CpAbsK=Up_fun(kk',meandepth,g,'');
CgAbsK=Ug_fun(kk,meandepth,g,'');

if strcmpi(input.wave.option,'No')
    if input.bdyassim.option==1
        nupeak=bdyassim.par.nupeak;
    else
        nupeak=ivp.par.nupeak;
    end
else
    nupeak=influx.par.nu_p;
end
nFreqP=model.nonlinear;

set(h.fig,'Renderer','zbuffer'); %due to graphics driver
%colorbar;
cstr={'*b','*g','*k','*c','*r'};
axis1=h.preview.dispersion_axes1;
hold(axis1,'off');

if strcmpi(model.dispersion(5:end),'Exact')
    plot(axis1,kk,OmegaAbsK,'r');
else
    OmExabsK=funOprt_OmExact(kk,meandepth,g,'');
    POm2=plot(axis1,kk,OmExabsK,'--r',kk,OmegaAbsK,'-r');
end

istr=1;
for ii=1:length(nupeak)
    if mod(istr,5)==0, istr=1;end
    for jj=1:nFreqP
        kpeak=invOm_fun(kk,jj*nupeak(ii),meandepth,g);
        hold(axis1,'on');
        plot(axis1,kpeak,jj*nupeak(ii),cstr{istr})
    end
    istr=istr+1;
end
xlim(axis1,[0 max(kk)]);
if ~strcmpi(model.dispersion(5:end),'Exact')
    legend(POm2,'Exact',model.dispersion(5:end),'location','NorthWest')
end
xlabel(axis1,'|k|');ylabel(axis1,'\omega [rad/s]')
grid(axis1,'off');
axes_properties(axis1,1.5);

ax_old = axis1;
f_new = figure('position',[1.5 1.5 0.7 0.7]);
set(f_new,'Renderer','zbuffer');
ax_new = copyobj(ax_old,f_new);
set(f_new,'visible','off')
axes(ax_new)
axes_properties(ax_new,1);
set(f_new,'unit','normalized','position',[1.5 1.5 0.7 0.7])
set(ax_new,'outerposition',[0 0 1 1]);
saveas(ax_new,[prevSaveDir,'dispersion.png']);
close(f_new)

axis2=h.preview.dispersion_axes2;
hold(axis2,'off');
if strcmpi(model.dispersion(5:end),'Exact')
    PC1=plot(axis2,kk,CpAbsK,'r',kk,CgAbsK,'b');
else
    CpExabsK=funOprt_UpExact(kk,meandepth,g);
    CgExabsK=funOprt_UgExact(kk,meandepth,g,'');
    PC2=plot(axis2,kk,CpAbsK,'r',kk,CpExabsK,'--r',kk,CgAbsK,'b',kk,CgExabsK,'--b');
end

istr=1;
for ii=1:length(nupeak)
    if mod(istr,5)==0, istr=1;end
    for jj=1:nFreqP
        hold(axis2,'on');
        kpeak=invOm_fun(kk,jj*nupeak(ii),meandepth,g);
        plot(axis2,kpeak,Up_fun(kpeak,meandepth,g,''),cstr{istr},...
            kpeak,Ug_fun(kpeak,meandepth,g,''),cstr{istr})
    end
    istr=istr+1;
end
hold(axis2,'off');
if strcmpi(model.dispersion(5:end),'Exact')
    legend(PC1,'Phase speed','Group speed');
else
    legend(PC2, ['Phase speed (',model.dispersion(5:end),')'],'Phase speed (Exact)',...
        ['Group speed (',model.dispersion(5:end),')'],'Group speed (Exact)')
end

xlim(axis2,[0 max(kk)]);
xlabel(axis2,'|k|');ylabel(axis2,'Speeds [m/s]')
grid(axis2,'off');
axes_properties(axis2,1.5);

ax_old = axis2;
f_new = figure('position',[1.5 1.5 0.7 0.7]);
set(f_new,'Renderer','zbuffer');
ax_new = copyobj(ax_old,f_new);
set(f_new,'visible','off')
axes(ax_new)
axes_properties(ax_new,1);
set(f_new,'unit','normalized','position',[1.5 1.5 0.7 0.7])
set(ax_new,'outerposition',[0 0 1 1]);
saveas(ax_new,[prevSaveDir,'speeds.png']);
close(f_new)



axes_properties(axis1,1.5);
axes_properties(axis2,1.5);
% for ii=1:8
%     if ii==8, id=5; else, id=ii; end
if strcmpi(input.wall.option,'Yes')
set(h.preview.spatial_popup_var,'value',5);
else
set(h.preview.spatial_popup_var,'value',1);
end
PreProcSpatialView(h,preproc);
% end
%%%%%%plot wave input
callback_preview_wave(h.preview.wave_popup_var,[],h);
%%%%%%%%%%%%%%%%%%%%%



%%%%log file
set(h.preview.log,'string',{log.string{:}})
end