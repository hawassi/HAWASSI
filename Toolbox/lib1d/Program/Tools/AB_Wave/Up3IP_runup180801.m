function [H_min,H_mid,H_plus,sg]= Up3IP_runup(k,Up,Om,omAdd,Cder,nupeak,cutfracwn,g,Depth,maxEtaInit,...
                            folder_root,savename,model)
% 3-parameter approx of phase velocity by interpolation
% cond: C, C0 and dC/deta exact exact at kapnu for disp of  C, C* and dC/deta operator
    
k_cut=max(k)/cutfracwn;
H_min  = (nupeak/k_cut)^2/g;%0.002;%
if model.nonlinear==2%strcmp(model.evol,'HSdirect')
H_plus = max(Depth(1),Depth(end))+maxEtaInit;    
H_mid  = (H_plus)*1/6;
elseif model.nonlinear>2 %strcmp(model.evol,'HS')
H_plus = max(Depth(1),Depth(end));    
H_mid  = (H_plus)*1/2;
end
Htot = linspace(H_min,H_plus,100);



% figure
% plot(Bathy)
% indxb
% indXs
if strcmp(model.dispersion,'OmSWE')
sg=[];
else
Ninterp=length(Htot);
Cp_H = zeros(Ninterp,1);
Cp_min = zeros(Ninterp,1);
Cp_plus = zeros(Ninterp,1);
Cp_mid = zeros(Ninterp,1);
C0_H = zeros(Ninterp,1);
C0_min = zeros(Ninterp,1);
C0_plus = zeros(Ninterp,1);
C0_mid = zeros(Ninterp,1);
Cder_H = zeros(Ninterp,1);
Cder_min = zeros(Ninterp,1);
Cder_plus = zeros(Ninterp,1);
Cder_mid = zeros(Ninterp,1);

k_half      =k(1:floor(end/2));
knu         = zeros(Ninterp,1);
for j= 1:Ninterp
Hj=Htot(j);
%kappa0_H= interp1(Om(k,Hj),k,k(1));
if strcmp(model.dispersion,'OmKdV')||strcmp(model.dispersion,'OmBBM')
    Omega       =Om(k_half,Hj,omAdd);
    [~,Idnu]    =max(Omega);
    kappanu_H   = interp1(Omega(1:Idnu),k_half(1:Idnu),nupeak,'spline');% invOmExact(2*pi/nupeak,depth);
else
    kappanu_H   = interp1(Om(k_half,Hj,omAdd),k_half,nupeak,'spline');% invOmExact(2*pi/nupeak,depth);
end
 knu(j)=  kappanu_H; 
%   if   nupeak<0.5
%       kappanu_H   = interp1(Om(k_half,Hj,omAdd),k_half,20*nupeak,'spline');%
%   end
%kappa_n_nu_H= k_max;%interp1(Om(k,Hj),k,nu);
%interpolation to find k at peak frequency nu

Cp_min(j) = Up(kappanu_H,H_min,omAdd);
Cp_plus(j) = Up(kappanu_H,H_plus,omAdd);
Cp_mid(j) = Up(kappanu_H,H_mid,omAdd);
Cp_H(j) = Up(kappanu_H,Hj,omAdd);

C0_min(j) = sqrt(g*H_min);
C0_plus(j)= sqrt(g*H_plus);
C0_mid(j) = sqrt(g*H_mid);
C0_H(j)   = sqrt(g*Hj);

Cder_min(j) = (Cder(kappanu_H,H_min));
Cder_plus(j) = (Cder(kappanu_H,H_plus));
Cder_mid(j) = (Cder(kappanu_H,H_mid));
Cder_H(j) = (Cder(kappanu_H,Hj));
end



A= (C0_mid.*Cder_plus)-(C0_plus.*Cder_mid);
B=-(C0_min.*Cder_plus)+(C0_plus.*Cder_min);
C= (C0_min.*Cder_mid)-(C0_mid.*Cder_min);
D=-(Cp_mid.*Cder_plus)+(Cp_plus.*Cder_mid);
E= (Cp_min.*Cder_plus)-(Cp_plus.*Cder_min);

F=-(Cp_min.*Cder_mid)+(Cp_mid.*Cder_min);
G= (Cp_mid.*C0_plus)-(Cp_plus.*C0_mid);
H=-(Cp_min.*C0_plus)+(Cp_plus.*C0_min);
I= (Cp_min.*C0_mid)-(Cp_mid.*C0_min);


det=Cp_min.*A+Cp_mid.*B+Cp_plus.*C;

gam_min  = ( A.*Cp_H+D.*C0_H+G.*Cder_H )./det;
gam_mid  = ( B.*Cp_H+E.*C0_H+H.*Cder_H )./det;
gam_plus = ( C.*Cp_H+F.*C0_H+I.*Cder_H )./det;



 if model.nonlinear==2
 sg.p=spline(Htot,gam_plus);
 sg.c=spline(Htot,gam_mid);
 sg.m=spline(Htot,gam_min);
 
 else
 sg.p=interp1(Htot,gam_plus,abs(Depth),'spline');
 sg.c=interp1(Htot,gam_mid,abs(Depth),'spline');
 sg.m=interp1(Htot,gam_min,abs(Depth),'spline');
end
 hfig=figure('Name','Bathy Set-up', 'Position',[170,10,2,600]);
 set(hfig,'Renderer','zbuffer'); %due to graphics driver
 subplot(4,1,1)
 plot(Htot,gam_min,'r',...
     Htot,gam_mid,'--g',...
     Htot,gam_plus,'-.b');
  title([savename,', ',model.dyn,', Interpol functions for minimal(b) and maximal depth vs depth']);
 %saveas(gcf,[folder_root,savename,dynmodel,'gamminmidplus.fig'])

 
%  pause

err   =zeros(length(Htot),length(k_half));
errknu=zeros(length(Htot),1);
Cnuv=zeros(length(Htot),1);
 for j=1:length(Htot)
     Cp_exact=Up(k_half,Htot(j),omAdd);
     Cp_approx=gam_min(j)*Up(k_half,H_min,omAdd)+...
                gam_plus(j)*Up(k_half,H_plus,omAdd)+...
                gam_mid(j)*Up(k_half,H_mid,omAdd); 
     err(j,:)=abs(Cp_exact-Cp_approx)./Cp_exact;
     indknu  =closest(k_half,knu(j));
     errknu(j)=abs(Cp_exact(indknu)-Cp_approx(indknu))/Cp_exact(indknu);
     Cnuv(j)=sqrt(Cp_exact(indknu));
 end
 
subplot(4,1,2)
mesh(Htot,k_half,err');
view(2);
xlabel('depth (m) ');
ylabel('wave number');
axis([0 max(Htot) 0 max(k_half)]);
hold on;
h=colorbar;
%caxis([0 1]);
title('Error of Interpolation of Cp')
ylabel(h, '|Exct-Approx|/Exct');
lineobj = findobj('type', 'line');
set(lineobj, 'linewidth', 1.3);
textobj = findobj('type', 'text');
set(textobj, 'fontunits', 'points');
set(textobj, 'fontsize', 10);
set(textobj, 'fontweight', 'bold');
%saveas(gcf,['ErrIn3p'],'fig');

 %Error at wavenumber in the peak frequency
subplot(4,1,3)
plot(Htot,errknu,'r');
title('Error of Cp (red) at peak frequencies')
xlabel('depth (m)');
ylabel('|Exct-Approx|/Exct');
grid on;    
lineobj = findobj('type', 'line');
set(lineobj, 'linewidth', 1.3);
textobj = findobj('type', 'text');
set(textobj, 'fontunits', 'points');
set(textobj, 'fontsize', 10);
set(textobj, 'fontweight', 'bold');
%saveas(gcf,['errknu'],'fig');          

subplot(4,1,4)
plot(Htot,knu);
title('wave number at peak frequency');
xlabel('depth');
ylabel('\kappa');


end


