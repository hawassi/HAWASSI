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

indxi=closest(x,Xinterv(1));
indxf =closest(x,Xinterv(end));
indti=closest(time,tinterv(1));
indtf =closest(time,tinterv(2));
timeNow=time(indti:dplot:indtf);
Nt=length(timeNow);
k=freqspace(x)';
g=9.81;
iter=1;Niter=Nt;tic;
set(jProgressBar,'Maximum',Niter, 'Value',iter);
statusbarObj.setText(['estimating time remaining..']);
Xinflux_R=Xinflux+0.1*lambda_p.*Nonlin_Adj;
Xinflux_L=Xinflux-0.1*lambda_p.*Nonlin_Adj;
indIL   =closest(x,Xinflux_R);
indFR   =closest(x,Xinflux_L);
indI     =indxi;
indF     =indxf;
MaxEtaInit=0.9*max(eta(:,closest(x,Xinflux)));
dx=x(2)-x(1);
depth=-bathy;
Ccrest=zeros(Nt,1);
Upart=zeros(Nt,1);
%figure
for i=1:Nt
    IDstop=eventLoopStop(jbStop);
    if IDstop==1, break;end;
    indt=indti+i*dplot-dplot;
    
    dx=x(2)-x(1);
    if any(eta(indt,[indI:indIL,indFR:indF])>MaxEtaInit);
        warning('off','all')
        [pks1,locs1]= findpeaks(eta(indt,[indI:indIL]),'minpeakheight',MaxEtaInit);
        [pks2,locs2]= findpeaks(eta(indt,[indFR:indF]),'minpeakheight',MaxEtaInit);
        indj1  =locs1(1:end)+indI-1;indj2  =locs2(1:end)+indFR-1;
        Nj1=length(indj1);Nj2=length(indj2);
        indj=zeros(Nj1+Nj2,1);pks=indj;
        indj(1:Nj1)=indj1;indj(Nj1+1:Nj1+Nj2)=indj2;
        pks(1:Nj1)=pks1;pks(Nj1+1:Nj1+Nj2)=pks2;
       % plot(x,eta(indt,:),'r',x([indj]),eta(indt,[indj]),'ob')

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
        Hx=imag(hilbert(eta(indt,:)));
        Kx=(eta(indt,:).*gradient(Hx,dx)-Hx.*gradient(eta(indt,:),dx))./(eta(indt,:).^2+Hx.^2);
        itCB=1;
        Totdepth=-bathy'+eta(indt,:);
       % for j=1:Npks
            [~,indjm]=max(pks);
            indjmax=indj(indjm);
            Kxj   =Kx(indjmax);
          
            Ccrest(iter)=abs(UpExact(Kxj,Totdepth(indjmax)));
            Upart(iter) =abs(u(indt,indjmax));
         %   pause
       % end
    end


if mod(iter,floor(0.1*Niter))==0 || iter==floor(0.01*Niter)
    set(jProgressBar,'Maximum',Niter, 'Value',iter);
    ETA=remain_time(iter,Niter);
    statusbarObj.setText(['time remaining=', num2str(ETA)]);
end
iter=iter+1;
end
Upart(Upart==0)=NaN;
Ccrest(Ccrest==0)=NaN;

set(jProgressBar,'Maximum',Niter, 'Value',Niter);
ETA=remain_time(Niter,Niter);
statusbarObj.setText(['time remaining=', num2str(ETA)]);
jProgressBar.setVisible(0);
jbStop.setVisible(0);
statusbarObj.setText('done.');

hf4=figure('Name','PostProc', 'Position',[250,10,800,500]); 
subplot(2,1,1)
plot(timeNow,Upart,'r',timeNow,Ccrest,'b');
legend('U','C')
xlabel('time');ylabel('Speeds');
title('Speeds at a maximum crest')
if GUIpp.PPsetting_xlim==1
    xlim([GUIpp.PPsetting_xlimVal(1) GUIpp.PPsetting_xlimVal(2)])
else
    xlim([timeNow(1) timeNow(iter-1)])
end
if GUIpp.PPsetting_ylim==1
    ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
end
plot_properties;

subplot(2,1,2)
plot(timeNow,Upart./Ccrest,'r');
xlabel('time');ylabel('U/C');
if GUIpp.PPsetting_xlim==1
    xlim([GUIpp.PPsetting_xlimVal(1) GUIpp.PPsetting_xlimVal(2)])
else
    xlim([timeNow(1) timeNow(iter-1)])
end
if GUIpp.PPsetting_ylim==1
    ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
end
plot_properties;
PPdata.PP1.Kinematic.time=timeNow;
PPdata.PP1.Kinematic.U=Upart;
PPdata.PP1.Kinematic.C=Ccrest;
PPdata.PP1.Kinematic.UperC=Upart./Ccrest;


