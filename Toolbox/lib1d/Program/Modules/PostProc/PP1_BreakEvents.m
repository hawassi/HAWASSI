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

tbreak =dataBreak_nodes(:,1);
iter=1;Niter=length(tbreak);tic;
set(jProgressBar,'Maximum',Niter, 'Value',iter);
statusbarObj.setText(['estimating time remaining..']);
xb=zeros(Niter,length(dataBreak_nodes(1,2:end)));
IdProp=sign(dataBreak_nodes(1,2)-dataBreak_nodes(1,3));
for I=1:length(tbreak)
    IDstop=eventLoopStop(jbStop);
    if IDstop==1, break;end;
    nodesBreak=[dataBreak_nodes(I,2:end) 0];
    [temp,locs]=findpeaks(nodesBreak);
    if IdProp==-1
    Nodes=[dataBreak_nodes(I,2) nodesBreak(locs+1)]; %take maxnodes+1 (maxnodes is trough)
    elseif IdProp==1 %(maxnodes)is crest
    Nodes=[dataBreak_nodes(I,2) temp];    
    end
    Nodes(Nodes==0)=1;
    xbreak=x(Nodes);
    xbreak(Nodes==1)=NaN;
    xb(I,1:length(xbreak))=xbreak;
    
    if mod(iter,floor(0.1*Niter))==0 || iter==floor(0.01*Niter)
        set(jProgressBar,'Maximum',Niter, 'Value',iter);
        ETA=remain_time(iter,Niter);
        statusbarObj.setText(['time remaining=', num2str(ETA)]);
    end
    iter=iter+1;
end
xb(xb==0)=NaN;
hf3=figure('Name','PostProc', 'Position',[250,10,800,500]);
plot(tbreak,xb,'ob');

ylabel('Breaking position [m]');
xlabel('time [s]');
PPdata.PP1.breaking.position=xb;
PPdata.PP1.breaking.time=tbreak;
plot_properties;

set(jProgressBar,'Maximum',Niter, 'Value',Niter);
ETA=remain_time(Niter,Niter);
statusbarObj.setText(['time remaining=', num2str(ETA)]);
jProgressBar.setVisible(0);
jbStop.setVisible(0);
statusbarObj.setText('done.');