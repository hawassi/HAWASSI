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

MWL_simul=mean(eta(indt1:dplot:indt2,indx1:indx2));
MTC_simul=max(eta(indt1:dplot:indt2,indx1:indx2));
MTT_simul=min(eta(indt1:dplot:indt2,indx1:indx2));
ATC_simul=zeros(size(MTC_simul));ATT_simul=ATC_simul;
iterx=1;
 timet=time(indt1:dplot:indt2);
for ii=indx1:indx2
    etat=eta(indt1:dplot:indt2,ii);
   
    datEtat=[timet etat];
    try
        [T, ind] = dat2wa(datEtat,MWL_simul(iterx),'c2c');
    catch
        continue;
    end
 etaTC=etat(ind);
 timeTC=timet(ind);
  timeC=timeTC(etaTC>MWL_simul(iterx));
  etaC=etaTC(etaTC>MWL_simul(iterx));
  etaT=etaTC(etaTC<MWL_simul(iterx));
  timeT=timeTC(etaTC<MWL_simul(iterx));
 
%  figure;
%  plot(timet,etat,'r',timeTC,etaTC,'or',timeC,etaC,'*b',timeT,etaT,'db')
%  
 ATC_simul(iterx)=mean(etaC);
 ATT_simul(iterx)=mean(etaT);
 iterx=iterx+1;
end
% 
% figure;
% plot(x(indx1:indx2),ATC_simul,'or',x(indx1:indx2),ATT_simul,'dr')

    if shippar.check==1
        for jj=1:shippar.Nship
            IndSL=closest(x,shippar.form.xShip(jj,1));
            IndSR=closest(x,shippar.form.xShip(jj,3));
            
            MTT_simul(IndSL+1:IndSR-1)=NaN;
            MTC_simul(IndSL+1:IndSR-1)=NaN; 
            ATC_simul(IndSL+1:IndSR-1)=NaN;
            ATT_simul(IndSL+1:IndSR-1)=NaN;
        end 
    end
          
          

hf=figure;
set(hf, 'units','normalized','Position', [0.2, 0.2, 0.4,0.4]);
subplot(2,1,1)
plot(x(indx1:indx2),MTC_simul,'-.k',x(indx1:indx2),MTT_simul,':c',...
    x(indx1:indx2),MWL_simul,'--r');
legend('MTC','MTT','MWL');
title('Maximal Temporal Amplitude');
xlabel('x[m]');ylabel('elevation [m]')
plot_properties;
subplot(2,1,2)
plot(x(indx1:indx2),ATC_simul,'-.k',x(indx1:indx2),ATT_simul,':c',...
    x(indx1:indx2),MWL_simul,'--r');
legend('ATC','ATT','MWL');
title('Average Temporal Amplitude');
if GUIpp.PPsetting_xlim==1
   ylim([GUIpp.PPsetting_xlimVal(1) GUIpp.PPsetting_xlimVal(2)]) 
else
xlim([x(indx1) x(indx2)])
end
if GUIpp.PPsetting_ylim==1
    ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
end
xlabel('x[m]');ylabel('elevation [m]')
plot_properties;



if GUIpp.save_figure==1
    if strcmp(GUIpp.save_figure_type,'.eps')
        saveas(hf,[sf_savename,'PP1_MTA_MWL',GUIpp.save_figure_type],'epsc')
    else
        saveas(hf,[sf_savename,'PP1_MTA_MWL',GUIpp.save_figure_type])
    end
end
