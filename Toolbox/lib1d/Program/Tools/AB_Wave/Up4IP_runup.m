function [sg]= Up4IP_runup(k,Up,Om,omAdd,nupeak,g,Depth,H_min,H_mid,H_plus,...
                            folder_root,savename,model)
% 3-parameter approx of phase velocity by interpolation
% cond: C, C0 and dC/deta exact exact at kapnu for disp of  C, C* and dC/deta operator
    

Htot = linspace(H_min,H_plus,100);
if strcmp(model.dispersion,'OmSWE')
sg=[];
else
Ninterp=length(Htot);
H_mid1=(H_min+H_mid)/2;
H_mid2=H_mid;%(H_mid+H_plus)/2;


k_half      =k(1:floor(end/2));
knu         = zeros(Ninterp,1);
gam_min     = zeros(Ninterp,1);
gam_mid1     = zeros(Ninterp,1);
gam_mid2     = zeros(Ninterp,1);
gam_plus    = zeros(Ninterp,1);

for j= 1:Ninterp
Hj=Htot(j);
%kappa0_H= interp1(Om(k,Hj),k,k(1));
if strcmp(model.dispersion,'OmKdV')||strcmp(model.dispersion,'OmBBM')
    Omega       =Om(k_half,Hj,omAdd);
    [~,Idnu]    =max(Omega);
    kappanu1_H   = interp1(Omega(1:Idnu),k_half(1:Idnu),nupeak,'spline');% invOmExact(2*pi/nupeak,depth);
else
    kappanuIG_H   = 0;%interp1(Om(k_half,Hj,omAdd),k_half,nupeak/10,'spline');
    kappanu1_H   = interp1(Om(k_half,Hj,omAdd),k_half,nupeak,'spline');% invOmExact(2*pi/nupeak,depth);
    kappanu2_H = interp1(Om(k_half,Hj,omAdd),k_half,nupeak*4,'spline');
    kappanu6_H = interp1(Om(k_half,Hj,omAdd),k_half,nupeak*6,'spline');
end
 knu(j)=  kappanu1_H;
%   if   nupeak<0.5
%       kappanu_H   = interp1(Om(k_half,Hj,omAdd),k_half,20*nupeak,'spline');%
%   end
%kappa_n_nu_H= k_max;%interp1(Om(k,Hj),k,nu);
%interpolation to find k at peak frequency nu

Cp1_min    = Up(kappanu1_H,H_min,omAdd).^2;
Cp1_plus   = Up(kappanu1_H,H_plus,omAdd).^2;
Cp1_mid1    = Up(kappanu1_H,H_mid1,omAdd).^2;
Cp1_mid2    = Up(kappanu1_H,H_mid2,omAdd).^2;
Cp1_H      = Up(kappanu1_H,Hj,omAdd).^2;

Cp2_min    = Up(kappanu2_H,H_min,omAdd).^2;
Cp2_plus   = Up(kappanu2_H,H_plus,omAdd).^2;
Cp2_mid1    = Up(kappanu2_H,H_mid1,omAdd).^2;
Cp2_mid2    = Up(kappanu2_H,H_mid2,omAdd).^2;
Cp2_H      = Up(kappanu2_H,Hj,omAdd).^2;

Cp6_min    = Up(kappanu6_H,H_min,omAdd).^2;% sqrt(g*H_min) ;%
Cp6_plus   = Up(kappanu6_H,H_plus,omAdd).^2;%%sqrt(g*H_plus) ;%
Cp6_mid1    = Up(kappanu6_H,H_mid1,omAdd).^2;%sqrt(g*H_mid) ;%
Cp6_mid2    = Up(kappanu6_H,H_mid2,omAdd).^2;%sqrt(g*H_mid) ;%
Cp6_H      = Up(kappanu6_H,Hj,omAdd).^2;% sqrt(g*Hj) ;%

CpIG_min    =  Up(kappanuIG_H,H_min,omAdd).^2;%
CpIG_plus   = Up(kappanuIG_H,H_plus,omAdd).^2;%%
CpIG_mid1    = Up(kappanuIG_H,H_mid1,omAdd).^2;%
CpIG_mid2    = Up(kappanuIG_H,H_mid2,omAdd).^2;%
CpIG_H      = Up(kappanuIG_H,Hj,omAdd).^2;%

C_Mat=[Cp1_min  Cp1_mid1  Cp1_mid2  Cp1_plus;...
       Cp2_min  Cp2_mid1  Cp2_mid2  Cp2_plus;...
       Cp6_min  Cp6_mid1  Cp6_mid2 Cp6_plus;...
       CpIG_min CpIG_mid1 CpIG_mid2 CpIG_plus];
C_vect=[Cp1_H; Cp2_H; Cp6_H; CpIG_H];
gam=C_Mat\C_vect;
gam_min(j)=gam(1);
gam_mid1(j)=gam(2);
gam_mid2(j)=gam(3);
gam_plus(j)=gam(4);
end


 if model.nonlinear==2
 sg.p=spline(Htot,gam_plus);
 sg.c1=spline(Htot,gam_mid1);
 sg.c2=spline(Htot,gam_mid2);
 sg.m=spline(Htot,gam_min);
 else
  DepthI=Depth;DepthI(Depth<H_min)=H_min;  
 sg.p=interp1(Htot,gam_plus,DepthI,'spline');
 sg.c=interp1(Htot,gam_mid1,DepthI,'spline');
 sg.m=interp1(Htot,gam_min,DepthI,'spline');
%  figure;
%  plot(Depth,sg.p,'r',Depth,sg.c,'g',Depth,sg.m,'b')
 end

 hfigInterp=figure('Name','Bathy Set-up', 'Position',[170,10,2,600]);
 set(hfigInterp,'Renderer','zbuffer'); %due to graphics driver
 subplot(4,1,1)
 plot([Htot],[gam_min.'],'r',...
     [Htot],[gam_mid1.'],'--g',...
     [Htot],[gam_mid2.'],'--b',...
     [Htot],[gam_plus.'],'-.k');
  title([savename,', ',model.dyn,', Interpol functions for minimal(b) and maximal depth vs depth']);
 %saveas(gcf,[folder_root,savename,dynmodel,'gamminmidplus.fig'])

 
%  pause
Htotnow=[linspace(0,H_min-H_min/10,20), Htot];
assignin('base','Htotnow',Htotnow)
err   =zeros(length(Htotnow),length(k_half));
errknu=zeros(length(Htotnow),1);

indHmin=closest(Htotnow,H_min);

knuj=zeros(size(Htotnow));
for j=1:length(Htotnow)
     Cp_exact=Up(k_half,Htotnow(j),omAdd).^2;
     if j<indHmin
     Cp_approx=(sqrt(g*Htotnow(j)).*ones(size(k_half))).^2;    
     else
     Cp_approx=gam_min(j-indHmin+1).*Up(k_half,H_min,omAdd).^2+...
                gam_plus(j-indHmin+1).*Up(k_half,H_plus,omAdd).^2+...
                gam_mid1(j-indHmin+1).*Up(k_half,H_mid1,omAdd).^2+...
                gam_mid2(j-indHmin+1).*Up(k_half,H_mid2,omAdd).^2; 
     end
     err(j,:)=abs(Cp_exact-Cp_approx)./Cp_exact;
    
     if j>1
     knuj(j)= interp1(Om(k_half,Htotnow(j),omAdd),k_half,nupeak,'spline');% 
     else
     knuj(j)=0.00001;    
     end
     indknu  =closest(k_half,knuj(j));
     errknu(j)=abs(Cp_exact(indknu)-Cp_approx(indknu))/Cp_exact(indknu);
 end
 
subplot(4,1,2)
mesh(Htotnow,k_half,err');
%mesh(Htotnow,k_half,errCder');
view(2);
xlabel('depth (m) ');
ylabel('wave number');
axis([0 max(Htotnow) 0 max(k_half)]);
hold on;
h=colorbar;
%caxis([0 0.5]);
title('Error of Cp^2 approx.')
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
plot(Htotnow,errknu,'r');
%plot(Htotnow,errknuCder,'r');
title('Error of Cp^2 (red) at peak wavenumbers')
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

% totDepthObs=[0.01;0.1;1;20;100];
% figure;
% colP={'r','g','b','c','y'};
% for ii=1:length(totDepthObs)
% wwj   = Om(k_half,totDepthObs(ii),omAdd);%
% Cpj   =Up(k_half,totDepthObs(ii),omAdd);
% plot(wwj,Cpj,colP{ii});
% hold on;
% end
% hold off;
% xlabel('\omega [rad/s]');
% ylabel('Cp');
% legend(['D=',num2str(totDepthObs(1))],['D=',num2str(totDepthObs(2))],...
%     ['D=',num2str(totDepthObs(3))],['D=',num2str(totDepthObs(4))],...
%     ['D=',num2str(totDepthObs(5))])
% plot_properties;
% 
% figure;
% for ii=1:length(totDepthObs)
% wwj   = Om(k_half,totDepthObs(ii),omAdd);%
% plot(k_half,wwj,colP{ii});
% hold on;
% end
% hold off;
% xlabel('k');
% ylabel('\omega [rad/s]');
% legend(['D=',num2str(totDepthObs(1))],['D=',num2str(totDepthObs(2))],...
%     ['D=',num2str(totDepthObs(3))],['D=',num2str(totDepthObs(4))],...
%     ['D=',num2str(totDepthObs(5))])
% plot_properties;
% 
% figure;
% for ii=1:length(totDepthObs)
% Cpj   = Up(k_half,totDepthObs(ii),omAdd);%
% plot(k_half,Cpj,colP{ii});
% hold on;
% end
% hold off;
% xlabel('k');
% ylabel('Cp ');
% legend(['D=',num2str(totDepthObs(1))],['D=',num2str(totDepthObs(2))],...
%     ['D=',num2str(totDepthObs(3))],['D=',num2str(totDepthObs(4))],...
%     ['D=',num2str(totDepthObs(5))])
% plot_properties;
