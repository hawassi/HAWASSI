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

if GUIpp.PPsetting_xcoarse==1
    dxplot= GUIpp.PPsetting_xcoarseVal;
else
    dxplot=1;
end

indxi=closest(x,Xinterv(1)); 
indxf =closest(x,Xinterv(end));
xNow=x(indxi:dxplot:indxf);
indti=closest(time,tinterv(1)); 
indtf =closest(time,tinterv(2));
timeNow=time(indti:dplot:indtf);
Nt=length(timeNow);
k=freqspace(x)';
g=9.81;
iter=1;Niter=Nt;tic;
set(jProgressBar,'Maximum',Niter, 'Value',iter);
statusbarObj.setText(['estimating time remaining..']);

BandFilt=GUIpp.Filter;
if BandFilt==1
    etaF=zeros(size(eta));
    uF=zeros(size(u));
    LF=GUIpp.LFHF(1);      %Low Frequency
    HF=GUIpp.LFHF(2);     %High Frequency
    statusbarObj.setText(['Filtering signals...']);
    for ii=1:length(xNow)
        indx=closest(x,xNow(ii));
        etaF(:,indx)=bandpass(time,eta(:,indx),LF,HF);
        uF(:,indx)=bandpass(time,u(:,indx),LF,HF);
    end
    eta=etaF;
    u=uF;
    statusbarObj.setText(['']);
end



if min(bathy)==max(bathy) %Flat bottom
  depth=-min(bathy);
  Cp=UpExact(k,depth)';
  Ep=zeros(Nt,1);
  Ek2=zeros(Nt,1);
  dx=x(2)-x(1);
if strcmp(dynmodel(1:3),'HS1')
    for i=1:Nt   
    IDstop=eventLoopStop(jbStop);
    if IDstop==1, break;end;
    
     indt=indti+i*dplot-dplot;
     u_hat=fft(u(indt,:));
     K2=Ifft(Cp.*u_hat).^2./g;
     Ek2(i)=trapz(K2(indxi:indxf))*dx/2;
     Ep(i)=trapz(g.*(eta(indt,indxi:indxf)).^2)*dx/2;
    
    if mod(iter,floor(0.1*Niter))==0 || iter==floor(0.01*Niter)
        set(jProgressBar,'Maximum',Niter, 'Value',iter);
        ETA=remain_time(iter,Niter);
        statusbarObj.setText(['time remaining=', num2str(ETA)]);
    end
    iter=iter+1;
    end
    Ek=Ek2;
    Etot=Ep+Ek;
    
     
elseif strcmp(dynmodel(1:3),'HS2')
    Ek3=zeros(Nt,1);
    M0=-1i.*k.*Cp.^2./g;
    for i=1:Nt
      IDstop=eventLoopStop(jbStop);
     if IDstop==1, break;end;   
     indt =indti+i*dplot-dplot;   
     u_now=u(indt,:);
     eta_now=eta(indt,:);
     u_hat=fft(u_now);
     Mu_hat=M0.*u_hat;
     Mu=Ifft(Mu_hat);
     K2=Ifft(Cp.*u_hat).^2./g;
     Ek2(i)=trapz(K2(indxi:indxf))*dx/2;
     K3=-Mu.^2.*eta_now+u_now.^2.*eta_now;
     Ek3(i)=trapz(K3(indxi:indxf))*dx/2;
     Ep(i)=trapz(g.*(eta(indt,indxi:indxf)).^2)*dx/2;
  
     if mod(iter,floor(0.1*Niter))==0 || iter==floor(0.01*Niter)
        set(jProgressBar,'Maximum',Niter, 'Value',iter);
        ETA=remain_time(iter,Niter);
        statusbarObj.setText(['time remaining=', num2str(ETA)]);
    end
    iter=iter+1;
    end
    Ek=Ek2+Ek3;
    Etot=Ep+Ek;
elseif strcmp(dynmodel(1:3),'HS3')
    Ek3=zeros(Nt,1);
    Ek4=zeros(Nt,1);
    M0=-1i.*k.*Cp.^2./g;
    for i=1:Nt
      IDstop=eventLoopStop(jbStop);
     if IDstop==1, break;end;    
     indt       =indti+i*dplot-dplot;
     u_now=u(indt,:);
     eta_now    =eta(indt,:);
     u_hat      =fft(u_now);
     Mu_hat     =M0.*u_hat;
     Mu         =Ifft(Mu_hat);
     dxetaMu_hat=1i.*k.*fft(Mu.*eta_now);
     dxetaMu    =Ifft(dxetaMu_hat);
     MdxetaMu_hat=M0.*dxetaMu_hat;
     MdxetaMu   =Ifft(MdxetaMu_hat);   
     dxu_hat    =1i.*k.*u_hat;
     dxu        =Ifft(dxu_hat);
     K2     =Ifft(Cp.*u_hat).^2./g;
     Ek2(i) =trapz(K2(indxi:indxf))*dx/2;
     K3     =-Mu.^2.*eta_now+u_now.^2.*eta_now;
     Ek3(i) =trapz(K3(indxi:indxf))*dx/2;
     K4     =eta_now.*Mu.*MdxetaMu+eta_now.^2.*dxu.*Mu;
     Ek4(i) =trapz(K4(indxi:indxf))*dx/2;
     Ep(i)=trapz(g.*eta_now(indxi:indxf).^2)*dx/2;
    
    
    if mod(iter,floor(0.1*Niter))==0 || iter==floor(0.01*Niter)
        set(jProgressBar,'Maximum',Niter, 'Value',iter);
        ETA=remain_time(iter,Niter);
        statusbarObj.setText(['time remaining=', num2str(ETA)]);
    end
    iter=iter+1;
    end
    Ek=Ek2+Ek3+Ek4;
    Etot=Ep+Ek;
elseif strcmp(dynmodel(1:3),'HS4')
    Ek3=zeros(Nt,1);
    Ek4=zeros(Nt,1);
    Ek5=zeros(Nt,1);
    M0=-1i.*k.*Cp.^2./g;
    for i=1:Nt
      IDstop=eventLoopStop(jbStop);
     if IDstop==1, break;end   
     indt       =indti+i*dplot-dplot;
     u_now      =u(indt,:);
     eta_now    =eta(indt,:);
     u_hat      =fft(u_now);
     Mu_hat     =M0.*u_hat;
     Mu         =Ifft(Mu_hat);
     dxMu_hat   =1i.*k.*Mu_hat;
     dxMu       =Ifft(dxMu_hat);
     dxxeta_hat =1i.*k.*1i.*k.*fft(eta_now);
     dxxeta     =Ifft(dxxeta_hat);
     dxetaMu_hat=1i.*k.*fft(Mu.*eta_now);
     dxetaMu    =Ifft(dxetaMu_hat);
     MdxetaMu_hat=M0.*dxetaMu_hat;
     MdxetaMu   =Ifft(MdxetaMu_hat);   
     dxu_hat    =1i.*k.*u_hat;
     dxu        =Ifft(dxu_hat);
     K2     =Ifft(Cp.*u_hat).^2./g;
     Ek2(i) =trapz(K2(indxi:indxf))*dx/2;
     K3     =-Mu.^2.*eta_now+u_now.^2.*eta_now;
     Ek3(i) =trapz(K3(indxi:indxf))*dx/2;
     K4     =eta_now.*Mu.*MdxetaMu+eta_now.^2.*dxu.*Mu;
     Ek4(i) =trapz(K4(indxi:indxf))*dx/2;
     K5     =-(1/3)*eta_now.^3.*(dxu).^2-eta_now.*MdxetaMu.^2 ...
             -eta_now.^2.*dxu.*MdxetaMu+(1/3).*eta_now.^3.*dxMu.^2 ...
             -(1/2).*eta_now.^2.*dxxeta.*Mu.^2;
     Ek5(i) =trapz(K5(indxi:indxf))*dx/2;    
     Ep(i)=trapz(g.*eta_now(indxi:indxf).^2)*dx/2;
   
    if mod(iter,floor(0.1*Niter))==0 || iter==floor(0.01*Niter)
        set(jProgressBar,'Maximum',Niter, 'Value',iter);
        ETA=remain_time(iter,Niter);
        statusbarObj.setText(['time remaining=', num2str(ETA)]);
    end
    iter=iter+1;
    end
    Ek=Ek2+Ek3+Ek4+Ek5;
    Etot=Ep+Ek;
end
elseif min(bathy)~=max(bathy) && ~any(bathy>0)%Varying bottom
   
  Ep=zeros(Nt,1);
  Ek2=zeros(Nt,1);
  dx=x(2)-x(1);
if strcmp(dynmodel(1:3),'HS1')
    for i=1:Nt   
    IDstop=eventLoopStop(jbStop);
    if IDstop==1, break;end;
    
     indt=indti+i*dplot-dplot;
     u_now=u(indt,:).';
     u_hat=fft(u_now);
     Cu_hat=OpCinterp(Oprt,u_hat);
     K2=Ifft(Cu_hat).^2./g;
     Ek2(i)=trapz(K2(indxi:indxf))*dx/2;
     Ep(i)=trapz(g.*(eta(indt,indxi:indxf)).^2)*dx/2;
    
    if mod(iter,floor(0.1*Niter))==0 || iter==floor(0.01*Niter)
        set(jProgressBar,'Maximum',Niter, 'Value',iter);
        ETA=remain_time(iter,Niter);
        statusbarObj.setText(['time remaining=', num2str(ETA)]);
    end
    iter=iter+1;
    end
    Ek=Ek2;
    Etot=Ep+Ek;
    
     
elseif strcmp(dynmodel(1:3),'HS2')
    Ek3=zeros(Nt,1);
     for i=1:Nt
      IDstop=eventLoopStop(jbStop);
     if IDstop==1, break;end;   
     indt =indti+i*dplot-dplot;   
     u_now=u(indt,:).';u_hat=fft(u_now);
     Cu_hat=OpCinterp(Oprt,u_hat);
     eta_now=eta(indt,:).';
     Mu_hat=Op_M0(k,g,Oprt,u_hat);
     Mu=Ifft(Mu_hat);
     K2=Ifft(Cu_hat).^2./g;
     Ek2(i)=trapz(K2(indxi:indxf))*dx/2;
     K3=-Mu.^2.*eta_now+u_now.^2.*eta_now;
     Ek3(i)=trapz(K3(indxi:indxf))*dx/2;
     Ep(i)=trapz(g.*(eta(indt,indxi:indxf)).^2)*dx/2;
  
     if mod(iter,floor(0.1*Niter))==0 || iter==floor(0.01*Niter)
        set(jProgressBar,'Maximum',Niter, 'Value',iter);
        ETA=remain_time(iter,Niter);
        statusbarObj.setText(['time remaining=', num2str(ETA)]);
    end
    iter=iter+1;
    end
    Ek=Ek2+Ek3;
    Etot=Ep+Ek;
elseif strcmp(dynmodel(1:3),'HS3')
    Ek3=zeros(Nt,1);
    Ek4=zeros(Nt,1);
    for i=1:Nt
      IDstop=eventLoopStop(jbStop);
     if IDstop==1, break;end    
     indt       =indti+i*dplot-dplot;
     u_now=u(indt,:).';u_hat      =fft(u_now);
     eta_now    =eta(indt,:).';
     Cu_hat=OpCinterp(Oprt,u_hat);
     Mu_hat     =Op_M0(k,g,Oprt,u_hat);
     Mu         =Ifft(Mu_hat);
     dxetaMu_hat=1i.*k.*fft(Mu.*eta_now);
     dxetaMu    =Ifft(dxetaMu_hat);
     MdxetaMu_hat=Op_M0(k,g,Oprt,dxetaMu_hat);
     MdxetaMu   =Ifft(MdxetaMu_hat);   
     dxu_hat    =1i.*k.*u_hat;
     dxu        =Ifft(dxu_hat);
     K2     =Ifft(Cu_hat).^2./g;
     Ek2(i) =trapz(K2(indxi:indxf))*dx/2;
     K3     =-Mu.^2.*eta_now+u_now.^2.*eta_now;
     Ek3(i) =trapz(K3(indxi:indxf))*dx/2;
     K4     =eta_now.*Mu.*MdxetaMu+eta_now.^2.*dxu.*Mu;
     Ek4(i) =trapz(K4(indxi:indxf))*dx/2;
     Ep(i)=trapz(g.*eta_now(indxi:indxf).^2)*dx/2;
    
    
    if mod(iter,floor(0.1*Niter))==0 || iter==floor(0.01*Niter)
        set(jProgressBar,'Maximum',Niter, 'Value',iter);
        ETA=remain_time(iter,Niter);
        statusbarObj.setText(['time remaining=', num2str(ETA)]);
    end
    iter=iter+1;
    end
    Ek=Ek2+Ek3+Ek4;
    Etot=Ep+Ek;
elseif strcmp(dynmodel(1:3),'HS4')
    Ek3=zeros(Nt,1);
    Ek4=zeros(Nt,1);
    Ek5=zeros(Nt,1);
    M0=-1i.*k.*Cp.^2./g;
    for i=1:Nt
      IDstop=eventLoopStop(jbStop);
     if IDstop==1, break;end;    
     indt       =indti+i*dplot-dplot;
     u_now      =u(indt,:).'; u_hat=fft(u_now);
     eta_now    =eta(indt,:).';
      Cu_hat=OpCinterp(Oprt,u_hat);
      Mu_hat     =Op_M0(k,g,Oprt,u_hat);
     Mu         =Ifft(Mu_hat);
     dxMu_hat   =1i.*k.*Mu_hat;
     dxMu       =Ifft(dxMu_hat);
     dxxeta_hat =1i.*k.*1i.*k.*fft(eta_now);
     dxxeta     =Ifft(dxxeta_hat);
     dxetaMu_hat=1i.*k.*fft(Mu.*eta_now);
     dxetaMu    =Ifft(dxetaMu_hat);
     MdxetaMu_hat=Op_M0(k,g,Oprt,dxetaMu_hat);
     MdxetaMu   =Ifft(MdxetaMu_hat);   
     dxu_hat    =1i.*k.*u_hat;
     dxu        =Ifft(dxu_hat);
     K2     =Ifft(Cu_hat).^2./g;
     Ek2(i) =trapz(K2(indxi:indxf))*dx/2;
     K3     =-Mu.^2.*eta_now+u_now.^2.*eta_now;
     Ek3(i) =trapz(K3(indxi:indxf))*dx/2;
     K4     =eta_now.*Mu.*MdxetaMu+eta_now.^2.*dxu.*Mu;
     Ek4(i) =trapz(K4(indxi:indxf))*dx/2;
     K5     =-(1/3)*eta_now.^3.*(dxu).^2-eta_now.*MdxetaMu.^2 ...
             -eta_now.^2.*dxu.*MdxetaMu+(1/3).*eta_now.^3.*dxMu.^2 ...
             -(1/2).*eta_now.^2.*dxxeta.*Mu.^2;
     Ek5(i) =trapz(K5(indxi:indxf))*dx/2;    
     Ep(i)=trapz(g.*eta_now(indxi:indxf).^2)*dx/2;
   
    if mod(iter,floor(0.1*Niter))==0 || iter==floor(0.01*Niter)
        set(jProgressBar,'Maximum',Niter, 'Value',iter);
        ETA=remain_time(iter,Niter);
        statusbarObj.setText(['time remaining=', num2str(ETA)]);
    end
    iter=iter+1;
    end
    Ek=Ek2+Ek3+Ek4+Ek5;
    Etot=Ep+Ek;
end

elseif any(bathy>0)%Run-Up

H_min      =OprtInterp_par.H_min1;
H_minShore  =OprtInterp_par.H_minShore;
IntCoef     =OprtInterp_par.IntCoef;

% [H_min,H_mid,H_plus,IntCoef]= Up3IP_runup(k,Up,Om,omAdd,Cder,nu_p,cutfrac,g,-bathy,...
%         maxEtaInit,[],[],model);
  

  Ep=zeros(Nt,1);
  Ek2=zeros(Nt,1);
  dx=x(2)-x(1);
  n  =length(x);


    Ek3=zeros(Nt,1);
  %  figure;
    for i=1:Nt
      IDstop=eventLoopStop(jbStop);
     if IDstop==1, break;end  
     indt =indti+i*dplot-dplot;
     
     eta_now    =eta(indt,:).';
     H=eta_now-bathy;
     
     HeavH=Heaviside(H-H_minShore);
     indDI=find(HeavH==0,1,'first');
     indDF=find(HeavH==0,1,'last');
     ShoreChar=HeavH;
     
     if isempty(indDI)
         indWI=1; indWF=n;
         indEndBreak=n;
         FlagOvertop=0;
         signProp=0;
     elseif indDI==1  % propagation direction to the left
         ShoreChar(indDI:indDF)=zeros;
         indWI=indDF+1; indWF=n;
         indEndBreak=indWI-10;
         FlagOvertop=0;
         signProp=-1;
     elseif indDF==n
         %     HeavHC2=Heaviside(H(indDI+1:indDF)-H_minShore);
         %     indDI=indDI+find(HeavHC2==0,1,'first');
         
         ShoreChar(indDI:indDF)=zeros;
         indWI=1; indWF=indDI-1;
         indEndBreak=indWF+10;
         FlagOvertop=0;
         signProp=1;
     elseif indDI~=1 && indDF~=n
         %HeavH(H<50*H_min)=0;
         indW=find(HeavH==1);
         indD=find(HeavH<H_minShore);
         ShoreChar=HeavH;
         indEndBreak=n;
         FlagOvertop=1;
         signProp=0;
     end

     if  strcmp(dynmodel(1:3),'HS2')
         C2m1         = Oprt.C2m1;
         C2p1         = Oprt.C2p1;
         C2c1         = Oprt.C2c1;
         C2m2         = Oprt.C2m2;
         C2p2         = Oprt.C2p2;
         C2c2         = Oprt.C2c2;
         
         gm1=zeros(n,1);gp1=gm1;gc1=gp1;
         gm2=zeros(n,1);gp2=gm2;gc2=gp2;
         gm1(indWI:indWF)=ppval(IntCoef.m1,H(indWI:indWF));
         gp1(indWI:indWF)=ppval(IntCoef.p1,H(indWI:indWF));
         gc1(indWI:indWF)=ppval(IntCoef.c1,H(indWI:indWF));
         gm2(indWI:indWF)=ppval(IntCoef.m2,H(indWI:indWF));
         gp2(indWI:indWF)=ppval(IntCoef.p2,H(indWI:indWF));
         gc2(indWI:indWF)=ppval(IntCoef.c2,H(indWI:indWF));

     else
%          Cm         = Up(k,H_min,omAdd)';
%          Cp         = Up(k,H_plus,omAdd)';
%          Cc         = Up(k,H_mid,omAdd)';
%      gm = IntCoef.m;   gp = IntCoef.p;  gc = IntCoef.c;  
     end
       
     u_now=u(indt,:).';u_now_hat=fft(u_now);
     
%      plot(x,u_now);
%      pause(0.001);
  
     Cu_hat= Op_C2r(sqrt(C2m1),gm1,sqrt(C2p1),gp1,sqrt(C2c1),gc1,...
                   sqrt(C2m2),gm2,sqrt(C2p2),gp2,sqrt(C2c2),gc2,u_now,u_now_hat);
      Cu=Ifft(Cu_hat);
      Cu(H<H_min) = sqrt(g*H(H<H_min)).*u_now(H<H_min);
               
     K2=Cu.^2./g;
   
     Ek2(i)=trapz(K2(indxi:indxf).*ShoreChar(indxi:indxf))*dx/2;
     Ep(i)=trapz(g.*(eta_now(indxi:indxf)).^2.*ShoreChar(indxi:indxf))*dx/2;
  
    
     
    if mod(iter,floor(0.1*Niter))==0 || iter==floor(0.01*Niter)
        set(jProgressBar,'Maximum',Niter, 'Value',iter);
        ETA=remain_time(iter,Niter);
        statusbarObj.setText(['time remaining=', num2str(ETA)]);
    end
    iter=iter+1;
    end
    Ek=Ek2+Ek3;
    Etot=Ep+Ek;


end
set(jProgressBar,'Maximum',Niter, 'Value',Niter); 
ETA=remain_time(Niter,Niter);
statusbarObj.setText(['time remaining=', num2str(ETA)]);
jProgressBar.setVisible(0);
jbStop.setVisible(0);
statusbarObj.setText('done.');

if IDstop==1, idtff=i-1;
else idtff=length(timeNow);end
hf2=figure('Name','PostProc', 'Position',[250,10,800,500]); 
plot(timeNow(1:idtff),Etot(1:idtff),'r',timeNow(1:idtff),Ep(1:idtff),'g',...
    timeNow(1:idtff),Ek(1:idtff),'--b')
legend('Total','Potential','Kinetic','Location','northwest');
xlabel('time[s]');
ylabel('Energy [J]');
if GUIpp.PPsetting_xlim==1
    xlim([GUIpp.PPsetting_xlimVal(1) GUIpp.PPsetting_xlimVal(2)])
else
    xlim([timeNow(1) timeNow(idtff)])
end
if GUIpp.PPsetting_ylim==1
    ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
end

plot_properties;
PPdata.PP1.Energy.Total=Etot;
PPdata.PP1.Energy.Potential=Ep;
PPdata.PP1.Energy.Potential=Ek;
