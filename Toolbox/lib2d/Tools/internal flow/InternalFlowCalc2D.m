function IFcalc=InternalFlowCalc2D(inputIF)
tic;
JavProgressBar;
[jbStop]=Java_stopbutton(statusbarObj);
jProgressBar.setStringPainted( true );
g=inputIF.data.par.g;
time=inputIF.data.outkinematic.time;
indtI=closest(time,inputIF.time.init);
indtF=closest(time,inputIF.time.end);
dtStep=inputIF.time.step;
T=time(indtI:dtStep:indtF);
indxI=closest(inputIF.data.dom.X,inputIF.x.init);
indxF=closest(inputIF.data.dom.X,inputIF.x.end);
dxStep=inputIF.x.step;
X=inputIF.data.dom.X(indxI:dxStep:indxF);
indyI=closest(inputIF.data.dom.Y,inputIF.y.init);
indyF=closest(inputIF.data.dom.Y,inputIF.y.end);
dyStep=inputIF.y.step;
Y=inputIF.data.dom.Y(indyI:dyStep:indyF);
dx=inputIF.data.dom.dx;
dy=inputIF.data.dom.dy;

if inputIF.z.Idequidistant==1
Z=linspace(inputIF.z.init,inputIF.z.end,inputIF.z.Nz);    
else
Z=inputIF.z.input;    
end
IFcalc.time=time(indtI:dtStep:indtF);
IFcalc.x=X;IFcalc.y=Y;IFcalc.z=Z;
bathy=inputIF.data.dom.bathy.profile;
Nx=length(X);Ny=length(Y);Nz=length(Z);Nt=length(T);
eta=inputIF.data.outkinematic.eta(indtI:dtStep:indtF,:,:);
phi=inputIF.data.outkinematic.phi(indtI:dtStep:indtF,:,:);
dteta=inputIF.data.outkinematic.dteta(indtI:dtStep:indtF,:,:);
dtphi=inputIF.data.outkinematic.dtphi(indtI:dtStep:indtF,:,:);

dtPhi=zeros(Nt,Nz,Ny,Nx);
dzPhi=zeros(Nt,Nz,Ny,Nx);
dyPhi=zeros(Nt,Nz,Ny,Nx);
dxPhi=zeros(Nt,Nz,Ny,Nx);
Ptot=zeros(Nt,Nz,Ny,Nx);
D=-bathy;
Hplus=max(max(max(eta)))+max(max(D));

if max(max(any(D<0)))==1
Hmin=inputIF.data.dom.bathy.HminShore;
IdShore=1;
else
Hmin=min(min(min(eta)))+min(min(D));
IdShore=0;
end
Htot=[Hmin:(Hplus-Hmin)/100:Hplus];
nupeak=inputIF.data.influx.par.nu_p;
KK=inputIF.data.dom.KK;
Kx=inputIF.data.dom.Kx;
Ky=inputIF.data.dom.Ky;
kk=sort(reshape(sqrt(Kx.^2+Ky.^2),[],1),'ascend');
absk=linspace(min(kk),max(kk),200).';

Omfun=str2func('OmExact');
[gam_Hmin0,gam_Hplus0]=sinhcosh_2IP(Htot,Hmin,Hplus,absk,Omfun,nupeak);
sgam_Hmin=spline(Htot,gam_Hmin0);
sgam_Hplus=spline(Htot,gam_Hplus0);

if min(min(D))==max(max(D))  %flat bottom
    ID_bathy='Flat';
else
    ID_bathy='NotFlat';
end

if ~strcmp(ID_bathy,'Flat')
    if IdShore==0
        Dmin=min(min(D));Dplus=max(max(D));
    else
        Dplus=max(max(D));Dmin=inputIF.data.dom.bathy.HminShore;
    end
    Dref=linspace(Dmin,Dplus,100);
    [gam_Dmin0,gam_Dplus0]=sinhcosh_2IP(Dref,Dmin,Dplus,absk,Omfun,nupeak);
    gam_Dmin=interp1(Dref,gam_Dmin0,D);
    gam_Dplus=interp1(Dref,gam_Dplus0,D);
    cosh_KD=gam_Dmin.*cosh(KK.*Dmin)+gam_Dplus.*cosh(KK.*Dplus);
    sinh_KD=gam_Dmin.*sinh(KK.*Dmin)+gam_Dplus.*sinh(KK.*Dplus);
end


iter=1;Niter=Nt;
set(jProgressBar,'Maximum',Niter, 'Value',iter);
statusbarObj.setText(['estimating time remaining ..']);


for i=1:Nt
    IDstop=eventLoopStop(jbStop);
    if IDstop==1, break;end
 phihat_i=fft2(squeeze(phi(i,:,:)));
 u       =funC_ifft2(1i.*Kx.*phihat_i);
 umax    =max(max(u));
 eta_i=squeeze(eta(i,:,:));

 dtphihat_i=fft2(squeeze(dtphi(i,:,:)));
 dteta_i=squeeze(dteta(i,:,:));
 
 
 H=(eta_i+D);

 gam_Hmin=ppval(sgam_Hmin,H);
 gam_Hplus=ppval(sgam_Hplus,H);

 cosh_KH=gam_Hmin.*cosh(KK.*Hmin)+gam_Hplus.*cosh(KK.*Hplus);
 sinh_KH=gam_Hmin.*sinh(KK.*Hmin)+gam_Hplus.*sinh(KK.*Hplus);
 tanh_KH=sinh_KH./cosh_KH;
  for l=1:Nz
     if strcmp(ID_bathy,'Flat') %&& Z(l)<min(Eta_i)% the condition z(l)<min(Eta_i) causes discontinuity for z>0 for pure flat bottom case
            cosh_KDZ=cosh(KK.*(Z(l)+D(1)));
            sinh_KDZ=sinh(KK.*(Z(l)+D(1)));
     else
       cosh_KDZ=cosh_KD.*cosh(KK.*Z(l))+sinh_KD.*sinh(KK.*Z(l));
       sinh_KDZ=sinh_KD.*cosh(KK.*Z(l))+cosh_KD.*sinh(KK.*Z(l));
     end
    
     COSH=cosh_KDZ./cosh_KH;
     SINH=sinh_KDZ./sinh_KH;
     factor=funC_ifft2(phihat_i.*COSH.*(-KK.*tanh_KH));
     
     tempdtPhi=funC_ifft2(dtphihat_i.*COSH)+dteta_i.*factor;
     tempdtPhi(Z(l)>eta_i)=NaN;
     tempdtPhi(Z(l)<bathy)=NaN;
     dtPhi(i,l,:,:)=tempdtPhi(indyI:dyStep:indyF,indxI:dxStep:indxF);
     
     tempdzPhi=funC_ifft2(phihat_i.*KK.*SINH);
     tempdzPhi(Z(l)>eta_i)=NaN;
     tempdzPhi(Z(l)<bathy)=NaN;
     dzPhi(i,l,:,:)=tempdzPhi(indyI:dyStep:indyF,indxI:dxStep:indxF);
     [dxH,dyH]=gradient(H,dx,dy);
     [dxD,dyD]=gradient(D,dx,dy);
     
     tempdxPhi=funC_ifft2(1i.*Kx.*phihat_i.*COSH)+dxH.*factor...
         +dxD.*tempdzPhi;
     tempdxPhi(Z(l)>eta_i)=NaN;
     tempdxPhi(Z(l)<bathy)=NaN;
    
     dxPhi(i,l,:,:)=tempdxPhi(indyI:dyStep:indyF,indxI:dxStep:indxF);
     
     tempdyPhi=funC_ifft2(1i.*Ky.*phihat_i.*COSH)+dyH.*factor...
         +dyD.*tempdzPhi;
     tempdyPhi(Z(l)>eta_i)=NaN;
     tempdyPhi(Z(l)<bathy)=NaN;
     dyPhi(i,l,:,:)=tempdyPhi(indyI:dyStep:indyF,indxI:dxStep:indxF);
     
     Ptot(i,l,:,:)=-dtPhi(i,l,:,:)-g.*Z(l)-0.5.*(dxPhi(i,l,:,:).^2+dyPhi(i,l,:,:).^2+dzPhi(i,l,:,:).*2); 
  end
  if mod(iter,floor(0.1*Niter))==0 || iter==floor(0.01*Niter)
      set(jProgressBar,'Maximum',Niter, 'Value',iter);
      ETA=remain_time(iter,Niter);
      statusbarObj.setText(['time remaining=', num2str(ETA)]);
  end
  iter=iter+1;
end

statusbarObj.setText(['Calculating acceleration...']);

[dt_dxPhi,dt_dyPhi,dt_dzPhi]=fun2D_fluid_acceleration(dxPhi,dyPhi,dzPhi,T,Z,Y,X);

IFcalc.dtPhi=dtPhi;
IFcalc.Ptot=Ptot;
IFcalc.dyPhi  =dyPhi;
IFcalc.dxPhi  =dxPhi;
IFcalc.dzPhi  =dzPhi;
%IFcalc.absVel=sqrt(dxPhi.^2+dyPhi.^2+dzPhi.^2);
IFcalc.dt_dxPhi=dt_dxPhi;
IFcalc.dt_dyPhi=dt_dyPhi;
IFcalc.dt_dzPhi=dt_dzPhi;
%IFcalc.abs_accel=sqrt(dt_dxPhi.^2+dt_dyPhi.^2+dt_dzPhi.^2);

Proj=inputIF.data.Proj;
statusbarObj.setText(['Saving data']);
save ('-v7.3',[Proj.workdir,Proj.savename,'_calc_Interior2D_',num2str(inputIF.fileId),'.mat'],'IFcalc')
     
set(jProgressBar,'Maximum',Niter, 'Value',Niter); 
ETA=remain_time(Niter,Niter);
statusbarObj.setText(['time remaining=', num2str(ETA)]);
jProgressBar.setVisible(0);
statusbarObj.setText('done.');
end