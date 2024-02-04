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

%%Kinematic Breaking Criteria
%%Wave will be break if particle speed exceeds crest speed
%%===============================================================

function [CB] = Kinematic_breaking_criterion(time,eta,u,x,depth,b,param,indEND)

global itCB
CB=[];
indI=param.xbstart;
if isempty(indEND)
indF=param.xbend;
else
indF=indEND;    
end

dx=x(2)-x(1);
if any(eta(indI:indF)>param.MaxEtaInit)
    warning('off','all')
    %        w = warning('query','last');  %turn off warning
    %        warning('off',w.identifier);
    if isempty(param.xbendL)&&isempty(param.xbstartR)
        try
            [pks,locs]= findpeaks(eta(indI:indF),'minpeakheight',param.MaxEtaInit,'MinPeakProminence',param.MaxEtaInit);
        catch
            [pks,locs]= findpeaks(eta(indI:indF),'minpeakheight',param.MaxEtaInit);
        end
        indj  =locs(1:end)+indI-1;
    else %in case of bi directional influx
        indIL=param.xbendL;indFR=param.xbstartR;
        try
            [pks1,locs1]= findpeaks(eta([indI:indIL]),'minpeakheight',param.MaxEtaInit,'MinPeakProminence',param.MaxEtaInit);
            [pks2,locs2]= findpeaks(eta([indFR:indF]),'minpeakheight',param.MaxEtaInit,'MinPeakProminence',param.MaxEtaInit);
        catch
            [pks1,locs1]= findpeaks(eta([indI:indIL]),'minpeakheight',param.MaxEtaInit);
            [pks2,locs2]= findpeaks(eta([indFR:indF]),'minpeakheight',param.MaxEtaInit);
        end
        
        indj1  =locs1(1:end)+indI-1;indj2  =locs2(1:end)+indFR-1;
        Nj1=length(indj1);Nj2=length(indj2);
        indj=zeros(Nj1+Nj2,1);pks=zeros(Nj1+Nj2,1);
        indj(1:Nj1)=indj1;indj(Nj1+1:Nj1+Nj2)=indj2;
        pks(1:Nj1)=pks1;pks(Nj1+1:Nj1+Nj2)=pks2;
    end
    
    if isempty(pks)
        flag=0;
    else
        flag=1;
    end
else
    flag=0;
end


if flag==1
    Npks=length(pks);
    Hx=imag(hilbert(eta));
    Kx=(eta.*gradient(Hx,dx)-Hx.*gradient(eta,dx))./(eta.^2+Hx.^2);
    itCB=1;
    
    Totdepth=depth+eta;
    for j=1:Npks
        indJ  =indj(j);
        Kxj   =Kx(indJ);
        Ccrest=UpExact(Kxj,Totdepth(indJ));
        Upart =u(indJ);
        
        if abs(Upart)>=(b*abs(Ccrest))
            CB.index(itCB)   =indJ;
            CB.Ucrest(itCB)  =abs(Upart);
            CB.DirProp(itCB) =sign(Upart);
            CB.position(itCB)=x(indJ);
            CB.time(itCB)    =time;
            itCB=itCB+1;
            %%        For checking
%                     if time>28
%                  %   figure
%                     subplot(4,1,3)
%                     plot(x,eta,'-k');
%                     xlim([20 100]);
%                     hold on;
%                     plot(x(indj),eta(indj),'ob');
%                     xlim([20 100]);
%                     title(['Profile,@ time',num2str(time)]);
%                     hold off;
%                     subplot(4,1,4)
%                     plot(x,u,'-k');
%                     xlim([20 100]);
%                     hold on;
%                     plot(x(indj),Ccrest,'or',x(indj),u(indj),'ob');
%                     title(['Horizontal Velocity,@ time',num2str(time)]);
%                     xlim([20 100]);
%                     hold off;
%                     end
        end
    end
end


