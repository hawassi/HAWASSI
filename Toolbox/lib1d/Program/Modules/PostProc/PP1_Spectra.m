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

JavProgressBar;
[jbStop]=Java_stopbutton(statusbarObj);
jProgressBar.setStringPainted( true );

Xinterv=GUIpp.PP1_Prof_X';
tinterv=GUIpp.PP1_Buoy_T;
if GUIpp.PPsetting_tcoarse==1
    dplot= GUIpp.PPsetting_tcoarseVal;
else
    dplot=1;
end

if GUIpp.PPsetting_xcoarse==1
    dxcoarse= GUIpp.PPsetting_xcoarseVal;
else
    dxcoarse=1;
end

indxi=closest(x,Xinterv(1));
indxf =closest(x,Xinterv(end));
indti=closest(time,tinterv(1));
indtf =closest(time,tinterv(2));
timeNow=time(indti:dplot:indtf);
Nt=length(timeNow);
xNow=x(indxi:dxcoarse:indxf);
Nx=length(xNow);

g=9.81;
iter=1;Niter=Nx;tic;
set(jProgressBar,'Maximum',Niter, 'Value',iter);
statusbarObj.setText(['estimating time remaining..']);

Spec_hat=zeros(floor(Nt/2),Nx);
    
%figure
for i=1:Nx
    IDstop=eventLoopStop(jbStop);
    if IDstop==1, break;end;
   
    simul=eta(indti:dplot:indtf,closest(x,xNow(i)));
   [OmS,Spec_hat(:,i)]=...
        funSP_variance_density_spectrum1d(timeNow,simul);
    
if mod(iter,floor(0.1*Niter))==0 || iter==floor(0.01*Niter)
    set(jProgressBar,'Maximum',Niter, 'Value',iter);
    ETA=remain_time(iter,Niter);
    statusbarObj.setText(['time remaining=', num2str(ETA)]);
end
iter=iter+1;
end

set(jProgressBar,'Maximum',Niter, 'Value',Niter);
ETA=remain_time(Niter,Niter);
statusbarObj.setText(['time remaining=', num2str(ETA)]);
jProgressBar.setVisible(0);
jbStop.setVisible(0);
statusbarObj.setText('done.');

[WW,XX]=meshgrid(OmS,xNow);

hf4=figure('Name','PostProc', 'Position',[250,10,800,500]); 

surf(XX,WW,Spec_hat.','edgecolor','none');
% hold on;
% [~,hc]=contour3(WW,XX,Spec_hat,'ShowText','on');
% set(hc,'linecolor','k','linewidth',2)    
%view(2);
xlabel('x [m]');ylabel('\omega [rad/s]');
zlabel('Var. dens. spect.')
if GUIpp.PPsetting_xlim==1
    xlim([GUIpp.PPsetting_xlimVal(1) GUIpp.PPsetting_xlimVal(2)])
    zlim2=max(max(Spec_hat(:,closest(xNow,GUIpp.PPsetting_xlimVal(1)):closest(xNow,GUIpp.PPsetting_xlimVal(2)))));
    zlim([0 zlim2]);caxis([0 zlim2]);
end
if GUIpp.PPsetting_ylim==1
    ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
    zlim2=max(max(Spec_hat(closest(OmS,GUIpp.PPsetting_ylimVal(1)):closest(OmS,GUIpp.PPsetting_ylimVal(2)),:)));
    zlim([0 zlim2]);caxis([0 zlim2]);
end
if GUIpp.PPsetting_xlim==1&& GUIpp.PPsetting_ylim==1
   zlim2=max(max(Spec_hat(closest(OmS,GUIpp.PPsetting_ylimVal(1)):closest(OmS,GUIpp.PPsetting_ylimVal(2)),closest(xNow,GUIpp.PPsetting_xlimVal(1)):closest(xNow,GUIpp.PPsetting_xlimVal(2)))));
   zlim([0 zlim2]); caxis([0 zlim2]);
end
shading interp;
plot_properties;



