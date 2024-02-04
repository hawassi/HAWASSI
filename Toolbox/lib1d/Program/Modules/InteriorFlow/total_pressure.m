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

function P_tot=total_pressure(IP)
dtPHI=IP.dtPHI;
Vx=IP.dxPHI;
Vz=IP.dzPHI;
X=IP.X;Z=IP.Z;
T=IP.timeIP;
P_tot=zeros(length(T),length(X),length(Z));
g=9.81;
for i=1:length(T)
 for l=1:length(Z) 
    P_tot(i,:,l)=-dtPHI(i,:,l)-g.*Z(l)-0.5.*(Vx(i,:,l).^2+Vz(i,:,l).*2);
 end   
end