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

if GUIpp.PPsetting_tcoarse==1
    dplot= GUIpp.PPsetting_tcoarseVal;
else
    dplot=1;
end

PPdata.PP1.shipRB=shipRB(1:dplot:end,:);
PPdata.PP1.ship_time=time(1:dplot:end);

Nship=shippar.Nship;

idCol=1;
CR={'r','--g','-.b',':m','or','sg','db','+m'};

hf5=figure('Name','PostProc: Ship displacement', 'Position',[250,10,800,500]); 

subplot(3,1,1)
for ii=1:Nship
 if idCol==9, idCol=1; end;
 plot(time,shipRB(:,ii),CR{idCol});
 hold on;
 idCol=idCol+1; 
 leg{ii}=['Ship #',num2str(ii)];
end
title(['Mean Surge:',num2str(mean(shipRB(:,1)),3),'std: ',num2str(std(shipRB(:,1)),3)])
xlabel('time [s]') ;   ylabel('X [m]');
if Nship>1
legend(leg{1:end},'Location','northwest');
end
tinterv=GUIpp.PP1_Buoy_T;
if GUIpp.PPsetting_xlim==1
    xlim([GUIpp.PPsetting_xlimVal(1) GUIpp.PPsetting_xlimVal(2)])
else
    xlim([tinterv(1) tinterv(2)])
end

if GUIpp.PPsetting_ylim==1
    ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
end
plot_properties;

subplot(3,1,2)
idCol=1;
for ii=1:Nship
 if idCol==9, idCol=1; end;
 plot(time,shipRB(:,2*Nship+ii),CR{idCol});
 hold on;
 idCol=idCol+1; 
 leg{ii}=['Ship #',num2str(ii)];
end
title('Heave')
xlabel('time [s]') ;   ylabel('Z [m]');
if Nship>1
legend(leg{1:end},'Location','northwest');
end
if GUIpp.PPsetting_xlim==1
    xlim([GUIpp.PPsetting_xlimVal(1) GUIpp.PPsetting_xlimVal(2)])
else
    xlim([tinterv(1) tinterv(2)])
end

if GUIpp.PPsetting_ylim==1
    ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
end
plot_properties;
assignin('base','shipRB',shipRB)
subplot(3,1,3)
idCol=1;
for ii=1:Nship
 if idCol==9, idCol=1; end
 plot(time,rad2deg(shipRB(:,4*Nship+ii)),CR{idCol});
 hold on;
 idCol=idCol+1; 
 leg{ii}=['Ship #',num2str(ii)];
end
title('Pitch')
xlabel('time [s]') ;   ylabel('Theta [deg]');
if Nship>1
legend(leg{1:end},'Location','northwest');
end
if GUIpp.PPsetting_xlim==1
    xlim([GUIpp.PPsetting_xlimVal(1) GUIpp.PPsetting_xlimVal(2)])
else
    xlim([tinterv(1) tinterv(2)])
end

if GUIpp.PPsetting_ylim==1
    ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
end
plot_properties;

smfact=1;%40;
[~,SpecX]=spect(time,shipRB(:,1),smfact,1);
[~,SpecZ]=spect(time,shipRB(:,3),smfact,1);
[OmS,SpecTheta]=spect(time,rad2deg(shipRB(:,5)),smfact,1);

%OmS=OmS/2/pi;%*10;

 SpecX=SpecX*2*pi;%/(10^5);
 SpecZ=SpecZ*2*pi;%/(10^5);
 SpecTheta=SpecTheta*2*pi;%/10;

figure;
subplot(3,1,1)
plot(OmS,SpecX,'r')
xlim([OmS(1) OmS(floor(end/2))]);
legend('Var. dens. sway')
ylabel('Resp. [m^2 s]')
%xlim([0 2.5])
plot_properties;
subplot(3,1,2)
plot(OmS,SpecZ,'r')
xlim([OmS(1) OmS(floor(end/2))])
ylabel('Resp. [m^2 s]')
legend('Var. dens. heave')
%xlim([0 2.5])

plot_properties;
subplot(3,1,3)
plot(OmS,SpecTheta,'r')
xlim([OmS(1) OmS(floor(end/2))])
ylabel('Resp. [deg^2 s]')
xlabel('\omega [rad/s]') ; 
legend('Var. dens. roll')
%xlim([0 2.5])
plot_properties;