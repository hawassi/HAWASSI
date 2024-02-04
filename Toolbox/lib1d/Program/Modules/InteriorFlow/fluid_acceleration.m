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

function [dt_dxPHI,dt_dzPHI]=fluid_acceleration(IP,iter,Niter)
%JavProgressBar;
%jProgressBar.setStringPainted( true ); 
%statusbarObj.setText(['estimating time remaining..']); 



dxPHI=IP.dxPHI;
dzPHI=IP.dzPHI;
Z    =IP.Z;
X    =IP.X;
T    =IP.timeIP;
dt_dxPHI=zeros(length(T),length(X),length(Z));
dt_dzPHI=zeros(length(T),length(X),length(Z));

for i=1:length(X)
    for j=1:length(Z)
     dt_dxPHI(2:end-1,i,j)=(dxPHI(3:end,i,j)-dxPHI(1:end-2,i,j))./(2*(T(2)-T(1)));
     dt_dzPHI(2:end-1,i,j)=(dzPHI(3:end,i,j)-dzPHI(1:end-2,i,j))./(2*(T(2)-T(1)));
    end
%     if mod(iter,floor(0.1*Niter))==0
%     set(jProgressBar,'Maximum',Niter, 'Value',iter);  
%     ETA=remain_time(iter,Niter);
%     statusbarObj.setText(['time remaining=', num2str(ETA)]);
%     end
% iter=iter+1;
end

% for i=1:length(T)
%     for j=1:length(X)
%      dt_dzPHI(i,j,2:end-1)=(dzPHI(i,j,3:end)-dzPHI(i,j,1:end-2))./(2*(T(2)-T(1)));
%     end
%     if mod(iter,floor(0.1*Ntot))==0
%     set(jProgressBar,'Maximum',Ntot, 'Value',iter);  
%     ETA=remain_time(iter,Ntot);
%     statusbarObj.setText(['time remaining=', num2str(ETA)]);
%     end
%     iter=iter+1;
% end

