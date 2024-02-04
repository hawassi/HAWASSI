function [H_min,sg,signProp]= Up2IP_runup(k,Up,Om,nupeak,cutfracwn,g,Bathy,maxEtaInit,...
                            folder_root,savename,dynmodel)
% 3-parameter approx of phase velocity by interpolation
% cond: C, C0 and dC/deta exact exact at kapnu for disp of  C, C* and dC/deta operator

k_cut=max(k)/cutfracwn;
H_min  = (nupeak/k_cut)^2/g;%0.002;%
H_plus = max(Bathy(1),Bathy(end))+maxEtaInit;

HeavH=Heaviside(Bathy-H_min);
if Bathy(end)<Bathy(1)
IndxH_min=find(HeavH==0,1,'first')-1;
indXs    = IndxH_min;  % index x at the shore
indxb=find(Bathy==Bathy(1),1,'last');
Ninterp=indXs-indxb+1;
indI     =indxb; indF=indXs;
else
IndxH_min=find(HeavH==0,1,'last')+1;
indXs    = IndxH_min;  % index x at the shore
indxb=find(Bathy==Bathy(end),1,'first');
Ninterp=-indXs+indxb+1; 
indI     =indXs; indF=indxb;
end
% figure
% plot(Bathy)
% indxb
% indXs


Htot=Bathy;
Cp_H = zeros(Ninterp,1);
Cp_min = zeros(Ninterp,1);
Cp_plus = zeros(Ninterp,1);
Cder_H = zeros(Ninterp,1);
Cder_min = zeros(Ninterp,1);
Cder_plus = zeros(Ninterp,1);

SIGN=sign(indXs-indxb);

for j= 1:Ninterp
Hj=Htot(indxb+(j-1).*SIGN);
%kappa0_H= interp1(Om(k,Hj),k,k(1));
kappanu_H= interp1(Om(k,Hj),k,nupeak);
if isnan(kappanu_H)
break; 
disp('adjust discretisation size');
end;
    
%kappa_n_nu_H= k_max;%interp1(Om(k,Hj),k,5*nu);

%interpolation to find k at peak frequency nu
Cp_min(j) = Up(kappanu_H,H_min);
Cp_plus(j) = Up(kappanu_H,H_plus);
Cp_H(j) = Up(kappanu_H,Hj);

Cder_min(j) = (Cder(kappanu_H,H_min));
Cder_plus(j) = (Cder(kappanu_H,H_plus));
Cder_H(j) = (Cder(kappanu_H,Hj));
end

    det         = Cp_min.*Cder_plus -Cp_plus.*Cder_min;
    gam_min     = ( Cder_plus.*Cp_H-Cp_plus.*Cder_H )./det;
    gam_plus    = ( -Cder_min.*Cp_H+Cp_min.*Cder_H)./det;



if SIGN==-1
gam_min=gam_min(Ninterp:-1:1);
gam_plus=gam_plus(Ninterp:-1:1);
end

 sg.p=spline(Htot(indI:indF),gam_plus);
 sg.m=spline(Htot(indI:indF),gam_min);
 
 signProp=SIGN;

figure
plot(Htot(indI:indF),gam_min,'r',...
     Htot(indI:indF),gam_plus,'-.b');
title([savename,', ',dynmodel,', Interpol functions for minimal(b) and maximal depth vs depth']);
%saveas(gcf,[folder_root,savename,dynmodel,'gamminmidplus.fig'])
