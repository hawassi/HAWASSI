function [sg]= Up2Derv3IP_runup(k,Up,Om,omAdd,detaC2,nupeak,g,Depth,H_min,H_mid,H_plus,...
                            folder_root,savename,model)
% 3-parameter approx of phase velocity by interpolation
% cond: C, C0 and dC/deta exact exact at kapnu for disp of  C, C* and dC/deta operator
    

Htot = linspace(H_min,H_plus,100);



% figure
% plot(Bathy)
% indxb
% indXs
if strcmp(model.dispersion,'OmSWE')
sg=[];
else
Ninterp=length(Htot);
detaC2k0_H = zeros(Ninterp,1);
detaC2k0_min = zeros(Ninterp,1);
detaC2k0_plus = zeros(Ninterp,1);
detaC2k0_mid = zeros(Ninterp,1);
detaC2kIG_H = zeros(Ninterp,1);
detaC2kIG_min = zeros(Ninterp,1);
detaC2kIG_plus = zeros(Ninterp,1);
detaC2kIG_mid = zeros(Ninterp,1);
detaC2kp_H = zeros(Ninterp,1);
detaC2kp_min = zeros(Ninterp,1);
detaC2kp_plus = zeros(Ninterp,1);
detaC2kp_mid = zeros(Ninterp,1);

k_half      =k(1:floor(end/2));
knu         = zeros(Ninterp,1);
% kappanu_H_min   = interp1(Om(k_half,H_min,omAdd),k_half,nupeak,'spline');%
% kappanu_H_plus   = interp1(Om(k_half,H_plus,omAdd),k_half,nupeak,'spline');%
% kappanu_H_mid   = interp1(Om(k_half,H_mid,omAdd),k_half,nupeak,'spline');%
% kappanuIG_H_min   = interp1(Om(k_half,H_min,omAdd),k_half,nupeak/4,'spline');%
% kappanuIG_H_plus   = interp1(Om(k_half,H_plus,omAdd),k_half,nupeak/4,'spline');%
% kappanuIG_H_mid   = interp1(Om(k_half,H_mid,omAdd),k_half,nupeak/4,'spline');%

for j= 1:Ninterp
Hj=Htot(j);
%kappa0_H= interp1(Om(k,Hj),k,k(1));
if strcmp(model.dispersion,'OmKdV')||strcmp(model.dispersion,'OmBBM')
    Omega       =Om(k_half,Hj,omAdd);
    [~,Idnu]    =max(Omega);
    kappanu_H   = interp1(Omega(1:Idnu),k_half(1:Idnu),nupeak,'spline');% invOmExact(2*pi/nupeak,depth);
else
    kappanu_H   = interp1(Om(k_half,Hj,omAdd),k_half,nupeak,'spline');% invOmExact(2*pi/nupeak,depth);
    kappanuIG_H   = interp1(Om(k_half,Hj,omAdd),k_half,nupeak/10,'spline');
    kappanu0_H = interp1(Om(k_half,Hj,omAdd),k_half,nupeak*3,'spline');
end
 knu(j)=  kappanu_H; 
%   if   nupeak<0.5
%       kappanu_H   = interp1(Om(k_half,Hj,omAdd),k_half,20*nupeak,'spline');%
%   end
%kappa_n_nu_H= k_max;%interp1(Om(k,Hj),k,nu);
%interpolation to find k at peak frequency nu

detaC2k0_min(j) = detaC2(kappanu0_H,H_min);
detaC2k0_plus(j)= detaC2(kappanu0_H,H_plus);
detaC2k0_mid(j) = detaC2(kappanu0_H,H_mid);
detaC2k0_H(j)   = detaC2(kappanu0_H,Hj);

detaC2kIG_min(j) = detaC2(kappanuIG_H,H_min);% sqrt(g*H_min) ;%
detaC2kIG_plus(j)= detaC2(kappanuIG_H,H_plus);%%sqrt(g*H_plus) ;%
detaC2kIG_mid(j) = detaC2(kappanuIG_H,H_mid);%sqrt(g*H_mid) ;%
detaC2kIG_H(j)   = detaC2(kappanuIG_H,Hj);% sqrt(g*Hj) ;%

detaC2kp_min(j) = (detaC2(kappanu_H,H_min));% Up(kappanuIG_H,H_min,omAdd);%
detaC2kp_plus(j)= (detaC2(kappanu_H,H_plus));%Up(kappanuIG_H,H_plus,omAdd);%%
detaC2kp_mid(j) = (detaC2(kappanu_H,H_mid));%Up(kappanuIG_H,H_mid,omAdd);%
detaC2kp_H(j)   = (detaC2(kappanu_H,Hj)); %Up(kappanuIG_H,Hj,omAdd);%
end


A= (detaC2kIG_mid.*detaC2kp_plus)-(detaC2kIG_plus.*detaC2kp_mid);
B=-(detaC2kIG_min.*detaC2kp_plus)+(detaC2kIG_plus.*detaC2kp_min);
C= (detaC2kIG_min.*detaC2kp_mid)-(detaC2kIG_mid.*detaC2kp_min);
D=-(detaC2k0_mid.*detaC2kp_plus)+(detaC2k0_plus.*detaC2kp_mid);
E= (detaC2k0_min.*detaC2kp_plus)-(detaC2k0_plus.*detaC2kp_min);

F=-(detaC2k0_min.*detaC2kp_mid)+(detaC2k0_mid.*detaC2kp_min);
G= (detaC2k0_mid.*detaC2kIG_plus)-(detaC2k0_plus.*detaC2kIG_mid);
H=-(detaC2k0_min.*detaC2kIG_plus)+(detaC2k0_plus.*detaC2kIG_min);
I= (detaC2k0_min.*detaC2kIG_mid)-(detaC2k0_mid.*detaC2kIG_min);


det=detaC2k0_min.*A+detaC2k0_mid.*B+detaC2k0_plus.*C;

gam_min  = ( A.*detaC2k0_H+D.*detaC2kIG_H+G.*detaC2kp_H )./det;
gam_mid  = ( B.*detaC2k0_H+E.*detaC2kIG_H+H.*detaC2kp_H )./det;
gam_plus = ( C.*detaC2k0_H+F.*detaC2kIG_H+I.*detaC2kp_H )./det;


 if model.nonlinear==2
 sg.p=spline(Htot,gam_plus);
 sg.c=spline(Htot,gam_mid);
 sg.m=spline(Htot,gam_min);
 
 else
  DepthI=Depth;DepthI(Depth<H_min)=H_min;  
 sg.p=interp1(Htot,gam_plus,DepthI,'spline');
 sg.c=interp1(Htot,gam_mid,DepthI,'spline');
 sg.m=interp1(Htot,gam_min,DepthI,'spline');
%  figure;
%  plot(Depth,sg.p,'r',Depth,sg.c,'g',Depth,sg.m,'b')
 end

 hfigInterp=figure('Name','Bathy Set-up', 'Position',[170,10,2,600]);
 set(hfigInterp,'Renderer','zbuffer'); %due to graphics driver
 subplot(4,1,1)
 plot(Htot,gam_min,'r',...
     Htot,gam_mid,'--g',...
     Htot,gam_plus,'-.b');
  title([savename,', ',model.dyn,', Interpol functions for minimal(b) and maximal depth vs depth']);
 %saveas(gcf,[folder_root,savename,dynmodel,'gamminmidplus.fig'])

 
%  pause
Htotnow=[linspace(0,H_min-H_min/10,20), Htot];
assignin('base','Htotnow',Htotnow)
errdetaC2   =zeros(length(Htotnow),length(k_half));
errknudetaC2=zeros(length(Htotnow),1);

indHmin=closest(Htotnow,H_min);
knuj=zeros(size(Htotnow));

detaC2_exactN=zeros(length(k_half),length(Htotnow));
detaC2_approxN=zeros(length(k_half),length(Htotnow));

for j=1:length(Htotnow)
      detaC2_exact=detaC2(k_half,Htotnow(j));
     if j<indHmin
     detaC2_approx=g.*ones(size(k_half));
     else
     detaC2_approx=gam_min(j-indHmin+1)*detaC2(k_half,H_min)+...
                gam_plus(j-indHmin+1)*detaC2(k_half,H_plus)+...
                gam_mid(j-indHmin+1)*detaC2(k_half,H_mid); 
     end
     errdetaC2(j,:)=abs(detaC2_exact-detaC2_approx)./detaC2_exact;
    
     if j>1
     knuj(j)= interp1(Om(k_half,Htotnow(j),omAdd),k_half,nupeak,'spline');% 
     else
     knuj(j)=0.00001;    
     end
     indknu  =closest(k_half,knuj(j));
     errknudetaC2(j)=abs(detaC2_exact(indknu)-detaC2_approx(indknu))/detaC2_exact(indknu);
      detaC2_exactN(:,j)=detaC2_exact.';
     detaC2_approxN(:,j)=detaC2_approx.';
 end
 
subplot(4,1,2)
%mesh(Htotnow,k_half,errdetaC2');
mesh(Htotnow,k_half,detaC2_exactN);
view(2);
xlabel('depth (m) ');
ylabel('wave number');
axis([0 max(Htot) 0 max(k_half)]);
hold on;
h=colorbar;
caxis([0 10]);
title('\partial_{\eta}C^2 exact.')
ylabel(h, '\partial_{\eta}C^2 exact.');
% title('Error of \partial_{\eta}C^2 approx.')
% ylabel(h, '|Exct-Approx|/Exct');
lineobj = findobj('type', 'line');
set(lineobj, 'linewidth', 1.3);
textobj = findobj('type', 'text');
set(textobj, 'fontunits', 'points');
set(textobj, 'fontsize', 10);
set(textobj, 'fontweight', 'bold');
%saveas(gcf,['ErrIn3p'],'fig');

 %Error at wavenumber in the peak frequency
subplot(4,1,3)
mesh(Htotnow,k_half,detaC2_approxN);
view(2);
xlabel('depth (m) ');
ylabel('wave number');
axis([0 max(Htotnow) 0 max(k_half)]);
hold on;
h=colorbar;
caxis([0 10]);
title('detaC2 approx.')
ylabel(h, 'detaC2 approx.');

% plot(Htotnow,errknudetaC2,'r');
% title('Error of \partial_{\eta}C^2 (red) at peak wavenumbers')
% xlabel('depth (m)');
% ylabel('|Exct-Approx|/Exct');
% grid on;    
lineobj = findobj('type', 'line');
set(lineobj, 'linewidth', 1.3);
textobj = findobj('type', 'text');
set(textobj, 'fontunits', 'points');
set(textobj, 'fontsize', 10);
set(textobj, 'fontweight', 'bold');
%saveas(gcf,['errknu'],'fig');          

subplot(4,1,4)
plot(Htotnow,knuj);
title('wave number at peak frequency');
xlabel('depth');
ylabel('\kappa');
saveas(hfigInterp,[folder_root,'Preview_interpOpBathy.fig'])
end

% Cp_H_app=gam_min.*Cp_min+gam_plus.*Cp_plus+gam_mid.*Cp_mid;
% C0_H_app=gam_min.*C0_min+gam_plus.*C0_plus+gam_mid.*C0_mid;
% Cder_H_app=gam_min.*Cder_min+gam_plus.*Cder_plus+gam_mid.*Cder_mid;

% figure;
% subplot(3,1,1)
% plot(Htot,Cp_min,'r',Htot,Cp_mid,'g',Htot,Cp_plus,'b',Htot,Cp_H,'k',Htot,Cp_H_app,'--c')
% xlabel('depth')
% legend('Cp_{min}','Cp_{mid}','Cp_{plus}','Cp_{Exact}','Cp_{Exact1}')
% subplot(3,1,2)
% plot(Htot,C0_min,'r',Htot,C0_mid,'g',Htot,C0_plus,'b',Htot,C0_H,'k',Htot,C0_H_app,'--c')
% xlabel('depth')
% legend('C0_{min}','C0_{mid}','C0_{plus}','C0_{Exact}','C0_{Exact1}')
% subplot(3,1,3)
% plot(Htot,Cder_min,'r',Htot,Cder_mid,'g',Htot,Cder_plus,'b',Htot,Cder_H,'k',Htot,Cder_H_app,'--c')
% xlabel('depth')
% legend('Cder_{min}','Cder_{mid}','Cder_{plus}','Cder_{Exact}','Cder_{Exact1}')
% 
% figure;
% subplot(2,1,1)
% plot(k_half,Up(k_half,H_min,omAdd),'r',k_half,Up(k_half,H_mid,omAdd),'g',k_half,Up(k_half,H_plus,omAdd),'b')
% legend('C(k,H_{min})','C(k,H_{mid})','C(k,H_{plus})')
% xlabel('k');ylabel('C')
% subplot(2,1,2)
% plot(k_half,Cder(k_half,H_min),'r',k_half,Cder(k_half,H_mid),'g',k_half,Cder(k_half,H_plus),'b')
% legend('\partial_{\eta}C(k,H_{min})','\partial_{\eta}C(k,H_{mid})','\partial_{\eta}C(k,H_{plus})')
% xlabel('k');ylabel('\partial_{\eta}C')

totDepthObs=[0.01;0.1;1;20;100];
figure;
colP={'r','g','b','c','y'};
for ii=1:length(totDepthObs)
wwj   = Om(k_half,totDepthObs(ii),omAdd);%
Cpj   =detaC2(k_half,totDepthObs(ii));
plot(wwj,Cpj,colP{ii});
hold on;
end
hold off;
xlabel('\omega [rad/s]');
ylabel('C\partial{\eta}C');
legend(['D=',num2str(totDepthObs(1))],['D=',num2str(totDepthObs(2))],...
    ['D=',num2str(totDepthObs(3))],['D=',num2str(totDepthObs(4))],...
    ['D=',num2str(totDepthObs(5))])
plot_properties;

figure;
for ii=1:length(totDepthObs)
wwj   = Om(k_half,totDepthObs(ii),omAdd);%
plot(k_half,wwj,colP{ii});
hold on;
end
hold off;
xlabel('k');
ylabel('\omega [rad/s]');
legend(['D=',num2str(totDepthObs(1))],['D=',num2str(totDepthObs(2))],...
    ['D=',num2str(totDepthObs(3))],['D=',num2str(totDepthObs(4))],...
    ['D=',num2str(totDepthObs(5))])
plot_properties;

figure;
for ii=1:length(totDepthObs)
Cpj   = detaC2(k_half,totDepthObs(ii));%
plot(k_half,Cpj,colP{ii});
hold on;
end
hold off;
xlabel('k');
ylabel('C\partial{\eta}C ');
legend(['D=',num2str(totDepthObs(1))],['D=',num2str(totDepthObs(2))],...
    ['D=',num2str(totDepthObs(3))],['D=',num2str(totDepthObs(4))],...
    ['D=',num2str(totDepthObs(5))])
plot_properties;
