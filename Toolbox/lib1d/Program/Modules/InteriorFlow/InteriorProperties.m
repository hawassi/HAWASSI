function [Ptot,dtPhi, dxPhi, dzPhi, dt_dxPHI, dt_dzPHI, Z, X,IDstop]=InteriorProperties(x,xinterv,dzOrZ,k,Om,nupeak,bathy,eta,u,dtphihat,dteta,timeIP,H_min,IDstop)
tic;
JavProgressBar;
[jbStop]=Java_stopbutton(statusbarObj);
jProgressBar.setStringPainted( true );

%[statusbarObj]=JavaFrame_handling();
g=9.81;
Nt=length(dtphihat(:,1));
D=-bathy;
sizD=size(D);
if sizD(1)>1
    D=D';
end
k=k';

Z=[dzOrZ(1):dzOrZ(2):dzOrZ(3)];

Nz=length(Z);
Indx1=closest(x,xinterv(1));
Indx2=closest(x,xinterv(2));
stepx=xinterv(3);
X=x(Indx1:stepx:Indx2);
dtPhi=zeros(Nt,length(X),Nz);
dzPhi=zeros(Nt,length(X),Nz);
dxPhi=zeros(Nt,length(X),Nz);
Ptot=zeros(Nt,length(X),Nz);
%dx_dtPHI=zeros(Nt,length(X),Nz);
dx=x(2)-x(1);
iter=1;Niter=Nt;
set(jProgressBar,'Maximum',Niter, 'Value',iter);
statusbarObj.setText(['estimating time remaining ..']);

if any(bathy>0) %runup case
    if (D(1))>0 %propagation to the right
        indxshore=find(D==0,1,'first');
        Hplus=max(max(eta(:,1:indxshore)))+max(D(1:indxshore));
    else %propagation to the left
        indxshore=find(D<0,1,'first')-1;
        Hplus=max(max(eta(:,indxshore:end)))+max(D(indxshore:end));
    end 
    Hmin=H_min;
else
    Hplus=max(max(eta))+max(D);
    Hmin=min(min(eta))+min(D);
end

Htot=[Hmin:(Hplus-Hmin)/100:Hplus];
[gam_Hmin0,gam_Hplus0]=sinhcosh_2IP(Htot,Hmin,Hplus,k,Om,nupeak);

sgam_Hmin=spline(Htot,gam_Hmin0);
sgam_Hplus=spline(Htot,gam_Hplus0);

if min(D)==max(D)  %flat bottom
    ID_bathy='Flat';
else
    ID_bathy='NotFlat';
    DZmin=0;DZplus=Hplus;
    DZtot=[DZmin:(DZplus-DZmin)/100:DZplus];
    assignin('base','DZtot',DZtot)
    [gam_DZmin0,gam_DZplus0]=sinhcosh_2IP(DZtot,DZmin,DZplus,k,Om,nupeak);
    sgam_DZmin=spline(DZtot,gam_DZmin0);
    sgam_DZplus=spline(DZtot,gam_DZplus0);
end

u_max=max(max(u));
assignin('base','u_max',u_max);
  %    fig11=   figure;
for i=1:Nt
    IDstop=eventLoopStop(jbStop);
    if IDstop==1, break;end;
    for l=1:Nz
        phihat=fft(u(i,:))./(1i.*k);
        H=(eta(i,:)+D);
        gam_Hmin=ppval(sgam_Hmin,H);
        gam_Hplus=ppval(sgam_Hplus,H);
        cosh_kH=gam_Hmin.*cosh(k.*Hmin)+gam_Hplus.*cosh(k.*Hplus);
        sinh_kH=gam_Hmin.*sinh(k.*Hmin)+gam_Hplus.*sinh(k.*Hplus);
        tanh_kH=sinh_kH./cosh_kH;

        if strcmp(ID_bathy,'Flat')
            cosh_kDZ=cosh(k.*(Z(l)+D(1)));
            sinh_kDZ=sinh(k.*(Z(l)+D(1)));
        else
            DZ=D+Z(l);
            DZplus=max(DZ);
            gam_DZmin=ppval(sgam_DZmin,DZ);
            gam_DZplus=ppval(sgam_DZplus,DZ);
            DZmin=min(DZ); 
            if DZmin<0, DZmin=0; end    
            cosh_kDZ=gam_DZmin.*cosh(k.*DZmin)+gam_DZplus.*cosh(k.*DZplus);
            sinh_kDZ=gam_DZmin.*sinh(k.*DZmin)+gam_DZplus.*sinh(k.*DZplus);
            
%             subplot(2,1,1)
%             plot(H,gam_Hmin,'or',H,gam_Hplus,'ob')
%             subplot(2,1,2)
%             plot(DZ,gam_DZmin,'or',DZ,gam_DZplus,'ob')
%             assignin('base','DZ',DZ)
%             pause(.001)
        end
        COSH=cosh_kDZ./cosh_kH;
        
        factor=Ifft(phihat.*COSH.*(-k.*tanh_kH));
        tempdtPhi=Ifft(dtphihat(i,:).*COSH)+dteta(i,:).*factor;
        %temp_dx_dtPHI=Ifft(1i.*k.*fft(tempdtPhi)); %%acceleration calculation using FFT
        tempdtPhi(Z(l)>eta(i,:))=NaN;
        tempdtPhi(Z(l)<bathy)=NaN;
        dtPhi(i,:,l)=tempdtPhi(Indx1:stepx:Indx2);
        tempdzPhi=Ifft(phihat.*k.*sinh_kDZ./cosh_kH);
        tempdzPhi(Z(l)>eta(i,:))=NaN;
        tempdzPhi(Z(l)<bathy)=NaN;
        dzPhi(i,:,l)=tempdzPhi(Indx1:stepx:Indx2);
        tempdxPhi=Ifft(1i.*k.*phihat.*COSH)+gradient(H,dx).*factor...
            +gradient(D,dx).*tempdzPhi;
        tempdxPhi(Z(l)>eta(i,:))=NaN;
        tempdxPhi(Z(l)<bathy)=NaN;
        dxPhi(i,:,l)=tempdxPhi(Indx1:stepx:Indx2);
        Ptot(i,:,l)=-dtPhi(i,:,l)-g.*Z(l)-0.5.*(dxPhi(i,:,l).^2+dzPhi(i,:,l).*2);
        %%acceleration calculation using FFT
        % temp_dx_dtPHI(Z(l)>eta(i,:))=NaN;
        % temp_dx_dtPHI(Z(l)<bathy)=NaN;
        % dx_dtPHI(i,:,l)=temp_dx_dtPHI(Indx1:stepx:Indx2);
    end
    
    if mod(iter,floor(0.1*Niter))==0 || iter==floor(0.01*Niter)
        set(jProgressBar,'Maximum',Niter, 'Value',iter);
        ETA=remain_time(iter,Niter);
        statusbarObj.setText(['time remaining=', num2str(ETA)]);
    end
    iter=iter+1;
end
IP.dxPHI  =dxPhi;
IP.dzPHI  =dzPhi;
IP.Z      =Z;
IP.X      =X;
IP.timeIP =timeIP;
if IDstop==0
statusbarObj.setText('please wait...');    
[dt_dxPHI, dt_dzPHI]=fluid_acceleration(IP,iter,Niter);
else
dt_dxPHI=[];
dt_dzPHI=[];
end

set(jProgressBar,'Maximum',Niter, 'Value',Niter); 
ETA=remain_time(Niter,Niter);
statusbarObj.setText(['time remaining=', num2str(ETA)]);
jProgressBar.setVisible(0);
statusbarObj.setText('done.');

