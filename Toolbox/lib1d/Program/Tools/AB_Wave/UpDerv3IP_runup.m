function [sg]= UpDerv3IP_runup(k,Up,Om,omAdd,Cder,nupeak,cutfracwn,g,Depth,H_min,H_mid,H_plus,...
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
Cderk0_H = zeros(Ninterp,1);
Cderk0_min = zeros(Ninterp,1);
Cderk0_plus = zeros(Ninterp,1);
Cderk0_mid = zeros(Ninterp,1);
CderkIG_H = zeros(Ninterp,1);
CderkIG_min = zeros(Ninterp,1);
CderkIG_plus = zeros(Ninterp,1);
CderkIG_mid = zeros(Ninterp,1);
Cderkp_H = zeros(Ninterp,1);
Cderkp_min = zeros(Ninterp,1);
Cderkp_plus = zeros(Ninterp,1);
Cderkp_mid = zeros(Ninterp,1);

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
    kappanu0_H = interp1(Om(k_half,Hj,omAdd),k_half,nupeak*4,'spline');
end
 knu(j)=  kappanu_H; 
%   if   nupeak<0.5
%       kappanu_H   = interp1(Om(k_half,Hj,omAdd),k_half,20*nupeak,'spline');%
%   end
%kappa_n_nu_H= k_max;%interp1(Om(k,Hj),k,nu);
%interpolation to find k at peak frequency nu

Cderk0_min(j) = Cder(kappanu0_H,H_min);
Cderk0_plus(j)= Cder(kappanu0_H,H_plus);
Cderk0_mid(j) = Cder(kappanu0_H,H_mid);
Cderk0_H(j)   = Cder(kappanu0_H,Hj);

CderkIG_min(j) = Cder(kappanuIG_H,H_min);% sqrt(g*H_min) ;%
CderkIG_plus(j)= Cder(kappanuIG_H,H_plus);%%sqrt(g*H_plus) ;%
CderkIG_mid(j) = Cder(kappanuIG_H,H_mid);%sqrt(g*H_mid) ;%
CderkIG_H(j)   = Cder(kappanuIG_H,Hj);% sqrt(g*Hj) ;%

Cderkp_min(j) = (Cder(kappanu_H,H_min));% Up(kappanuIG_H,H_min,omAdd);%
Cderkp_plus(j)= (Cder(kappanu_H,H_plus));%Up(kappanuIG_H,H_plus,omAdd);%%
Cderkp_mid(j) = (Cder(kappanu_H,H_mid));%Up(kappanuIG_H,H_mid,omAdd);%
Cderkp_H(j)   = (Cder(kappanu_H,Hj)); %Up(kappanuIG_H,Hj,omAdd);%
end


A= (CderkIG_mid.*Cderkp_plus)-(CderkIG_plus.*Cderkp_mid);
B=-(CderkIG_min.*Cderkp_plus)+(CderkIG_plus.*Cderkp_min);
C= (CderkIG_min.*Cderkp_mid)-(CderkIG_mid.*Cderkp_min);
D=-(Cderk0_mid.*Cderkp_plus)+(Cderk0_plus.*Cderkp_mid);
E= (Cderk0_min.*Cderkp_plus)-(Cderk0_plus.*Cderkp_min);

F=-(Cderk0_min.*Cderkp_mid)+(Cderk0_mid.*Cderkp_min);
G= (Cderk0_mid.*CderkIG_plus)-(Cderk0_plus.*CderkIG_mid);
H=-(Cderk0_min.*CderkIG_plus)+(Cderk0_plus.*CderkIG_min);
I= (Cderk0_min.*CderkIG_mid)-(Cderk0_mid.*CderkIG_min);


det=Cderk0_min.*A+Cderk0_mid.*B+Cderk0_plus.*C;

gam_min  = ( A.*Cderk0_H+D.*CderkIG_H+G.*Cderkp_H )./det;
gam_mid  = ( B.*Cderk0_H+E.*CderkIG_H+H.*Cderkp_H )./det;
gam_plus = ( C.*Cderk0_H+F.*CderkIG_H+I.*Cderkp_H )./det;


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
errCder   =zeros(length(Htotnow),length(k_half));
errknuCder=zeros(length(Htotnow),1);

Cder_exactN=zeros(length(k_half),length(Htotnow));
Cder_approxN=zeros(length(k_half),length(Htotnow));


indHmin=closest(Htotnow,H_min);
knuj=zeros(size(Htotnow));
for j=1:length(Htotnow)
      Cder_exact=Cder(k_half,Htotnow(j));
     if j<indHmin
     Cder_approx=(g./2)./sqrt(g.*Htotnow(j)).*ones(size(k_half));
     else
     Cder_approx=gam_min(j-indHmin+1)*Cder(k_half,H_min)+...
                gam_plus(j-indHmin+1)*Cder(k_half,H_plus)+...
                gam_mid(j-indHmin+1)*Cder(k_half,H_mid); 
     end
     errCder(j,:)=abs(Cder_exact-Cder_approx);
    
     if j>1
     knuj(j)= interp1(Om(k_half,Htotnow(j),omAdd),k_half,nupeak,'spline');% 
     else
     knuj(j)=0.00001;    
     end
     indknu  =closest(k_half,knuj(j));
     errknuCder(j)=abs(Cder_exact(indknu)-Cder_approx(indknu));
     Cder_exactN(:,j)=Cder_exact.';
     Cder_approxN(:,j)=Cder_approx.';
 end
 
subplot(4,1,2)
mesh(Htotnow,k_half,errCder');
mesh(Htotnow,k_half,Cder_exactN);
view(2);
xlabel('depth (m) ');
ylabel('wave number');
axis([0 max(Htot) 0 max(k_half)]);
hold on;
h=colorbar;
caxis([0 2]);
title(' \partial_{\eta}C exact')
ylabel(h, '\partial_{\eta}C exact.');
lineobj = findobj('type', 'line');
set(lineobj, 'linewidth', 1.3);
textobj = findobj('type', 'text');
set(textobj, 'fontunits', 'points');
set(textobj, 'fontsize', 10);
set(textobj, 'fontweight', 'bold');
%saveas(gcf,['ErrIn3p'],'fig');

 %Error at wavenumber in the peak frequency
subplot(4,1,3)
 mesh(Htotnow,k_half,Cder_approxN);
view(2);
xlabel('depth (m) ');
ylabel('wave number');
axis([0 max(Htot) 0 max(k_half)]);
hold on;
h=colorbar;
caxis([0 2]);
title(' \partial_{\eta}C approx')
ylabel(h, '\partial_{\eta}C approx.');

% plot(Htotnow,errknuCder,'r');
% title('Error of \partial_{\eta}C (red) at peak wavenumbers')
% xlabel('depth (m)');
% ylabel('|Exct-Approx|');
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

% totDepthObs=[0.01;0.1;1;20;100];
% figure;
% colP={'r','g','b','c','y'};
% for ii=1:length(totDepthObs)
% wwj   = Om(k_half,totDepthObs(ii),omAdd);%
% Cpj   =Cder(k_half,totDepthObs(ii));
% plot(wwj,Cpj,colP{ii});
% hold on;
% end
% hold off;
% xlabel('\omega [rad/s]');
% ylabel('\partial{\eta}C');
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
% Cpj   = Cder(k_half,totDepthObs(ii));%
% plot(k_half,Cpj,colP{ii});
% hold on;
% end
% hold off;
% xlabel('k');
% ylabel('\partial{\eta}C ');
% legend(['D=',num2str(totDepthObs(1))],['D=',num2str(totDepthObs(2))],...
%     ['D=',num2str(totDepthObs(3))],['D=',num2str(totDepthObs(4))],...
%     ['D=',num2str(totDepthObs(5))])
% plot_properties;
