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
tsnap=GUIpp.PP1_Prof_T;

if GUIpp.PPsetting_xcoarse==1
    dplot= GUIpp.PPsetting_xcoarseVal;
else
    dplot=1;
end


if strcmp(wall.presence,'Yes')
    Xwall=wall.position; 
    ampli=max(max(eta));
end

dt=time(2)-time(1);

% if  strcmp(dynmodel,'HS1FWall')||strcmp(dynmodel,'HS2FWall')||strcmp(dynmodel,'HS3FWall')||strcmp(dynmodel,'HS4FWall')
%     if isempty(handles.input)  || Load_data_Id==1;
%         par=GUIpp.MySim.par;
%         influx_dat=GUIpp.MySim.influx;
%     else
%         par=evalin('base','par');
%         influx_dat=evalin('base','influx');
%     end
%     if isfield(par,'wall')
%     reflPercent=par.wall.refl_Coef*100;
%     Xwall=par.wall.position;    
%     else
%     reflPercent=par.reflPercent;
%     Xwall=par.wallposition;
%     end
%     
%     kp=influx_dat.k_p;
%     ampli=influx_dat.Hs./2;
%     TT=influx_dat.Tp;
%     Xinflux=influx_dat.position;
%     lambda_nu=2*pi/kp;
%     Xwall0=x(closest(x,Xwall));
%     Xinflux=x(closest(x,Xinflux)-1);
%     tsnap=time(closest(time,tsnap));
%     dt=time(2)-time(1);
%     hold on;
%     tsnapW=tsnap
%     x_inp=x;dx=x(2)-x(1);
%     theta1= kp.*(x_inp-Xinflux)-tsnapW.*2*pi./TT;
%     theta2=-kp.*(x_inp-(Xinflux+2*(Xwall0-Xinflux)))-tsnapW.*2*pi./TT;
% %         if strcmp(dynmodel,'HS1FWall')
%         Eta_Sol1=ampli.*cos(theta1);
%         Eta_Sol2=(reflPercent/100).*ampli.*cos(theta2);
% 
% %         else
% %         depth=-bathy(1);
% %         sigma=tanh(kp*depth);
% %         Eta_Sol1=ampli.*(cos(theta1)+kp.*ampli.*(3-sigma^2)./(4*sigma^3)*cos(2*theta1));%+ampli.*cos(2*Xwall+Xinflux-xinp-t_inp.*2*pi/TT);
% %         Eta_Sol2=(reflPercent/100).*ampli.*(cos(theta2)+kp.*(reflPercent/100).*ampli.*(3-sigma^2)./(4*sigma^3)*cos(2*theta1));%+ampli.*cos(2*Xwall+Xinflux-xinp-t_inp.*2*pi/TT)
% %         end
% %     
%     indXwall=closest(x,Xwall);
%     indXinflux=closest(x,Xinflux);
%     Eta_Sol1(indXwall+1:end)=zeros;
%     Eta_Sol2(indXwall+1:end)=zeros;
%     Eta_Sol1(1:indXinflux-1)=zeros;
%     Eta_Sol2(1:indXinflux-1)=zeros;
%     Eta_Sol=Eta_Sol1+Eta_Sol2;
%     plot(x,Eta_Sol,'--b')
%     
% end


if ~isempty(dataBreak_nodes)
   PP1_Profile_breaking;
else    
   PP1_Profile_nonbreaking;
end
plot_properties;

xlabel('x [m]') ; ylabel('\eta [m]');
if GUIpp.PPsetting_xlim==1
    xlim([GUIpp.PPsetting_xlimVal(1) GUIpp.PPsetting_xlimVal(2)])
else
    xlim([Xinterv(1) Xinterv(2)])
end
if GUIpp.PPsetting_ylim==1
    ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
end

grid on;