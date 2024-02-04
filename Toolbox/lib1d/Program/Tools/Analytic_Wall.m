%%Analytic Wall
clc;close all;
xx=0:0.1:160;
depth=5;
ampl=0.05;Tp=2;
ww=2*pi/Tp;
kk=invOmExact(ww,depth);
Xwall=149.9997;
dx=xx(2)-xx(1);

Ywall=linspace(-2*ampl,2*ampl,100);
Xwall=xx(closest(xx,Xwall));
xinf=0.0244;
xinf1=xinf;
xinf2=Xwall+Xwall-xinf;
fig=figure;
%%%wave-profile
tinterv=0:0.1:600;
tsnap=522.2;%200:0.1:600;
for ii=1:length(tsnap)
eta_prof_L=ampl*(cos(kk.*(xx-xinf1)-tsnap(ii)*ww));
eta_prof_L(xx>Xwall)=0;
eta_prof_R=ampl*(cos(-kk.*(xx-xinf2)-tsnap(ii)*ww));
eta_prof_LR=eta_prof_L+eta_prof_R;
eta_prof_LR(xx>Xwall)=0;
plot(xx,eta_prof_L,'r',xx,eta_prof_R,'--b',xx,eta_prof_LR,'-.k',Xwall.*ones(size(Ywall)),Ywall,'or');
title(['t= ',num2str(tsnap(ii)),'[s]'])
xlim([100 150])
pause(0.1)
end