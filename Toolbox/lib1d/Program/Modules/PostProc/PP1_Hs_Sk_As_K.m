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

indx1=closest(x,Xinterv(1));
indx2 =closest(x,Xinterv(2));
indt1=closest(time,tinterv(1));
indt2 =closest(time,tinterv(2));
Asym_simul=zeros(1,length(indx1:indx2));
Skew_simul=Asym_simul;Hs_simul=Asym_simul;Kurt_simul=Hs_simul;

iter=1;Niter=indx2-indx1;tic;
set(jProgressBar,'Maximum',Niter, 'Value',iter);
statusbarObj.setText(['estimating time remaining..']);


for i=indx1:indx2
     IDstop=eventLoopStop(jbStop);
    if IDstop==1, break;end;
    
    MeanSimul=mean(eta(indt1:dplot:indt2,i));
    Simul=eta(indt1:indt2,i);
    Hs_simul(i-indx1+1)= 4*sqrt(var(Simul-MeanSimul,1));
    Hilb_sim=imag(hilbert(Simul-MeanSimul));
    Asym_simul(i-indx1+1)=mean(Hilb_sim.^3)/((mean((Simul-MeanSimul).^2))^(3/2));
    Skew_simul(i-indx1+1)=mean((Simul-MeanSimul).^3)/((mean((Simul-MeanSimul).^2))^(3/2));
    Kurt_simul(i-indx1+1)=kurtosis(Simul);
    
    if mod(iter,floor(0.1*Niter))==0 || iter==floor(0.01*Niter)
        set(jProgressBar,'Maximum',Niter, 'Value',iter);
        ETA=remain_time(iter,Niter);
        statusbarObj.setText(['time remaining=', num2str(ETA)]);
    end
    iter=iter+1;
end


if GUIpp.PPsetting_xcoarse==1
    stepdx= GUIpp.PPsetting_xcoarseVal;
else
    stepdx=1;
end

hf=figure;
set(hf, 'units','normalized','Position', [0.2, 0.2, 0.4,0.4]);
subplot(4,1,1)
plot(x(indx1:stepdx:indx2),Hs_simul(1:stepdx:end),'r')
xlabel('x[m]');
ylabel('Hs [m]')
plot_properties;
xlim([x(indx1) x(indx2)])
if GUIpp.PPsetting_ylim==1
    ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
end
subplot(4,1,2)
plot(x(indx1:stepdx:indx2),Skew_simul(1:stepdx:end),'r')
xlabel('x[m]');
ylabel('Sk ')
plot_properties;
xlim([x(indx1) x(indx2)])
if GUIpp.PPsetting_ylim==1
    ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
end
subplot(4,1,3)
plot(x(indx1:stepdx:indx2),-Asym_simul(1:stepdx:end),'r')
xlabel('x[m]');
ylabel('-Asym')
plot_properties;
xlim([x(indx1) x(indx2)])
if GUIpp.PPsetting_ylim==1
    ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
end
subplot(4,1,4)
plot(x(indx1:stepdx:indx2),Kurt_simul(1:stepdx:end),'r')
xlabel('x[m]');
ylabel('Kurt')
plot_properties;
xlim([x(indx1) x(indx2)])
if GUIpp.PPsetting_ylim==1
    ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
end


if GUIpp.save_figure==1
    if strcmp(GUIpp.save_figure_type,'.eps')
        saveas(hf,[sf_savename,'PP3_Hs_Sk_Asym',GUIpp.save_figure_type],'epsc')
    else
        saveas(hf,[sf_savename,'PP3_Hs_Sk_Asym',GUIpp.save_figure_type])
    end
end

set(jProgressBar,'Maximum',Niter, 'Value',Niter); 
ETA=remain_time(Niter,Niter);
statusbarObj.setText(['time remaining=', num2str(ETA)]);
jProgressBar.setVisible(0);
jbStop.setVisible(0);
statusbarObj.setText('done.');