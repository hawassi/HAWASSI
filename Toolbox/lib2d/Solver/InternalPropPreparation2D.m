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

%%%%%%%%%%%%%%%%Internal Flow module Preparation%%%%%%%%%%%%%%%%%%%%%%%%%%%
global   dteta dtphi iterInterior timeIP phiIP etaIP
dteta=[]; dtphi=[]; iterInterior=1; timeIP=[]; phiIP=[]; etaIP=[];

if options.interior.check==1
    
    
    if tintval(1)> par.IP.time(2) || tintval(end)< par.IP.time(1)
        checkInterior=0;
    else
        checkInterior=1;
        if Npartition==1
            t_IP_start=par.IP.time(1);
            t_IP_end  =par.IP.time(2);
        else
            if tintval(1)> par.IP.time(1)
                t_IP_start=tintval(1);
            else
                t_IP_start=par.IP.time(1);
            end
            
            if tintval(end)< par.IP.time(2)
                t_IP_end=tintval(end);
            else
                t_IP_end=par.IP.time(2);
            end
            
        end
        
        
        Nt_IP=round((t_IP_end-t_IP_start)/(par.IP.time(3).*timeSimul.dt));
        dteta=zeros(Nt_IP,dom.Ny,dom.Nx);
        dtphi=zeros(Nt_IP,dom.Ny,dom.Nx);
        timeIP  =zeros(Nt_IP,1);
        phiIP=zeros(Nt_IP,dom.Ny,dom.Nx);
        etaIP=zeros(Nt_IP,dom.Ny,dom.Nx);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
