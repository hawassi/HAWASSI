function [D_min,D_plus,gam_min,gam_plus]= Up2IP(k,Up,Om,omAdd,nupeak,Depth,...
                            folder_root,savename,dynmodel)
% 2-parameter approx of phase velocity by interpolation
% (from Ruddy Beji-Battjes): cond: Om exact and C^2 exact at kapnu for disp of  C^2
% operator
g=9.81;
nu      = nupeak;
D_min   = min(Depth);
D_plus  = max(Depth);
DDepth   =linspace(D_min,D_plus,100);

kappanu_D   = zeros(length(DDepth),1);
Cp_D        = zeros(length(DDepth),1);
Cp_min      = zeros(length(DDepth),1);
Cp_plus     = zeros(length(DDepth),1);
Om_D        = zeros(length(DDepth),1);
Om_min      = zeros(length(DDepth),1);
Om_plus     = zeros(length(DDepth),1);
% assignin('base','DDepth',DDepth)
for j= 1:length(DDepth)
    D               = DDepth(j);
    kappanu_D(j)    = funOprt_invOmExactVal(nu,D,g);%interp1(Om(k,D,omAdd),k,nu);%interpolation to find k at peak frequency nu
    Cp_min(j)       = Up(kappanu_D(j),D_min,omAdd).^2;
    Cp_plus(j)      = Up(kappanu_D(j),D_plus,omAdd).^2;
    Cp_D(j)         = Up(kappanu_D(j),D,omAdd).^2;
    Om_min(j)       = Om(kappanu_D(j),D_min,omAdd);
    Om_plus(j)      = Om(kappanu_D(j),D_plus,omAdd);
    Om_D(j)         = Om(kappanu_D(j),D,omAdd);
end
% assignin('base','Om_D',Om_D);
% assignin('base','Om_min',Om_min);
% assignin('base','Om_plus',Om_plus);


    det         = Cp_min.*Om_plus -Cp_plus.*Om_min;
    Ggam_min     = ( Om_plus.*Cp_D-Cp_plus.*Om_D )./det;
    Ggam_plus    = ( -Om_min.*Cp_D+Cp_min.*Om_D)./det;
    
 gam_plus=interp1(DDepth,Ggam_plus,Depth,'spline');
 gam_min=interp1(DDepth,Ggam_min,Depth,'spline');
       
 
% figure(105)
% plot(Bathy,gam_min,'b', Bathy, gam_plus,'r')  
% title([savename,', ',dynmodel,', Interpol functions for minimal(b) and maximal depth vs depth']);
% saveas(gcf,[folder_root,savename,dynmodel,'gamminplus.fig'])
%%

%%
% end
hfig=figure('Name','Bathy Set-up', 'Position',[170,10,2,600]);
 set(hfig,'Renderer','zbuffer'); %due to graphics driver
 subplot(4,1,1)
 plot(DDepth,Ggam_min,'r',...
     DDepth,Ggam_plus,'-.b');
title([savename,', ',dynmodel,', Interpol functions for minimal(b) and maximal depth vs depth']);
 %saveas(gcf,[folder_root,savename,dynmodel,'gamminmidplus.fig'])

%  pause
k_half=k(1:floor(end/16));
err   =zeros(length(DDepth),length(k_half));
errknu=zeros(length(DDepth),1);
Cnuv=zeros(length(DDepth),1);
 for j=1:length(DDepth)
     Cp_exact=Up(k_half,DDepth(j),omAdd).^2;
     Cp_approx=Ggam_min(j)*Up(k_half,D_min,omAdd).^2+...
                Ggam_plus(j)*Up(k_half,D_plus,omAdd).^2; 
     err(j,:)=abs(Cp_exact-Cp_approx)./abs(Cp_exact);
     indknu  =closest(k_half,kappanu_D(j));
     errknu(j)=abs(Cp_exact(indknu)-Cp_approx(indknu))/abs(Cp_exact(indknu));
     Cnuv(j)=sqrt(Cp_exact(indknu));
 end
 
subplot(4,1,2)
mesh(DDepth,k_half,err');
view(2);
xlabel('depth (m) ');
ylabel('wave number');
axis([min(DDepth) max(DDepth) 0 max(k_half)/2]);
hold on;
h=colorbar;
%caxis([0 1]);
%title('Error of Interpolation of Cp')
ylabel(h, '|Exact-Approx|/|Exact|');
lineobj = findobj('type', 'line');
set(lineobj, 'linewidth', 1.3);
textobj = findobj('type', 'text');
set(textobj, 'fontunits', 'points');
set(textobj, 'fontsize', 12);
set(textobj, 'fontweight', 'bold');
%saveas(gcf,['ErrIn3p'],'fig');
plot_properties;
set(textobj, 'fontsize', 18);
set(lineobj, 'linewidth', 2.3);
 %Error at wavenumber in the peak frequency
subplot(4,1,3)
plot(DDepth,errknu,'r');
%title('Error of Cp (red) at peak frequencies')
xlabel('depth (m)');
ylabel('|Exact-Approx|/|Exact|');
xlim([min(DDepth) max(DDepth)]);
grid on;    
lineobj = findobj('type', 'line');
set(lineobj, 'linewidth', 1.3);
textobj = findobj('type', 'text');
set(textobj, 'fontunits', 'points');
set(textobj, 'fontsize', 12);
set(textobj, 'fontweight', 'bold');
plot_properties;

set(textobj, 'fontsize', 18);
set(lineobj, 'linewidth', 2.3);
%saveas(gcf,['errknu'],'fig');          

subplot(4,1,4)
plot(DDepth,kappanu_D);
title('wave number at peak frequency');
xlabel('depth');
ylabel('\kappa');