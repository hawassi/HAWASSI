function [InterpVg]= funOprt_GroupVelInterpolation_3p(g,dom,Ug_fun,Om_fun,omAdd,nupeak,Depth,bath,...
                            Proj,dynmodel)
% 2-parameter approx of dispersion relation by interpolation
% cond: Cp  and Om^2 exact at kapnu for d Om^2 operator

nu      = mean(nupeak);
D_min   = min(min(Depth));
D_plus  = max(max(Depth));
D_mid   = (D_plus+D_min)/2;
DDepth   =linspace(D_min,D_plus,100);

kappanu_D   = zeros(length(DDepth),1);
Cg_D        = zeros(length(DDepth),1);
Cg_min      = zeros(length(DDepth),1);
Cg_mid      = zeros(length(DDepth),1);
Cg_plus     = zeros(length(DDepth),1);
Om_D        = zeros(length(DDepth),1);
Om_min      = zeros(length(DDepth),1);
Om_mid      = zeros(length(DDepth),1);
Om_plus     = zeros(length(DDepth),1);
Cp_D        = zeros(length(DDepth),1);
Cp_min      = zeros(length(DDepth),1);
Cp_mid      = zeros(length(DDepth),1);
Cp_plus     = zeros(length(DDepth),1);


if strcmp(bath.type,'Bx')
k=dom.kx;
else
k=dom.ky;    
end

for j= 1:length(DDepth)
    D               = DDepth(j);
    kappanu_D(j)    = funOprt_invOmExactVal(nu,D,g);%interp1(Om_fun(k,D,g),k,nu);%interpolation to find k at peak frequency nu
    Cg_min(j)       = Ug_fun(kappanu_D(j),D_min,g,omAdd);
    Cg_mid(j)       = Ug_fun(kappanu_D(j),D_mid,g,omAdd);
    Cg_plus(j)      = Ug_fun(kappanu_D(j),D_plus,g,omAdd);
    Cg_D(j)         = Ug_fun(kappanu_D(j),D,g,omAdd);
    Om_min(j)       = Om_fun(kappanu_D(j),D_min,g,omAdd);
    Om_mid(j)       = Om_fun(kappanu_D(j),D_mid,g,omAdd);
    Om_plus(j)      = Om_fun(kappanu_D(j),D_plus,g,omAdd);
    Om_D(j)         = Om_fun(kappanu_D(j),D,g,omAdd);
    Cp_min(j)       = funOprt_UpExact(kappanu_D(j),D_min,g).^2;
    Cp_mid(j)       = funOprt_UpExact(kappanu_D(j),D_mid,g).^2;
    Cp_plus(j)      = funOprt_UpExact(kappanu_D(j),D_plus,g).^2;
    Cp_D(j)         = funOprt_UpExact(kappanu_D(j),D,g).^2;
end

A= (Om_mid.*Cp_plus)-(Om_plus.*Cp_mid);
B=-(Om_min.*Cp_plus)+(Om_plus.*Cp_min);
C= (Om_min.*Cp_mid)-(Om_mid.*Cp_min);
D=-(Cg_mid.*Cp_plus)+(Cg_plus.*Cp_mid);
E= (Cg_min.*Cp_plus)-(Cg_plus.*Cp_min);

F=-(Cg_min.*Cp_mid)+(Cg_mid.*Cp_min);
G= (Cg_mid.*Om_plus)-(Cg_plus.*Om_mid);
H=-(Cg_min.*Om_plus)+(Cg_plus.*Om_min);
I= (Cg_min.*Om_mid)-(Cg_mid.*Om_min);


det=Cg_min.*A+Cg_mid.*B+Cg_plus.*C;


Ggam_min  = ( A.*Cg_D+D.*Om_D+G.*Cp_D )./det;
Ggam_mid  = ( B.*Cg_D+E.*Om_D+H.*Cp_D )./det;
Ggam_plus = ( C.*Cg_D+F.*Om_D+I.*Cp_D )./det;

 gam_plus=interp1(DDepth,Ggam_plus,Depth,'spline');
 gam_min=interp1(DDepth,Ggam_min,Depth,'spline');
 gam_mid=interp1(DDepth,Ggam_mid,Depth,'spline');
 InterpVg.gam_plus=gam_plus;
 InterpVg.gam_min= gam_min;
 InterpVg.gam_mid= gam_mid;
 InterpVg.D_min  =D_min;
 InterpVg.D_mid  =D_mid;
 InterpVg.D_plus =D_plus;
 InterpVg.DDepth =DDepth;
 InterpVg.Ggam_min =Ggam_min;
 InterpVg.Ggam_plus =Ggam_plus;
 InterpVg.Ggam_mid =Ggam_mid;

 
%  [X,Y]=meshgrid(dom.X,dom.Y);
%  figure;
%  surf(X,Y,gam_min,'edgecolor','none')
       
% figure(105)
% plot(Bathy,gam_min,'b', Bathy, gam_plus,'r')  
% title([savename,', ',dynmodel,', Interpol functions for minimal(b) and maximal depth vs depth']);
% saveas(gcf,[folder_root,savename,dynmodel,'gamminplus.fig'])
%%

%%
% end
%  figure
%  subplot(4,1,1)
%  plot(DDepth,Ggam_min,'r',DDepth,Ggam_mid,'-.g',...
%      DDepth,Ggam_plus,':b');
% % title([Proj.savename,', ',dynmodel,', Interpol functions for minimal(b) and maximal depth vs depth']);
%  %saveas(gcf,[folder_root,savename,dynmodel,'gamminmidplus.fig'])
% 
% %  pause
% k_half=k(1:floor(end/8))';
% err   =zeros(length(DDepth),length(k_half));
% errknu=zeros(length(DDepth),1);
% Cnuv=zeros(length(DDepth),1);
%  for j=1:length(DDepth)
%      Cp_exact=Ug_fun(k_half,DDepth(j),omAdd);
%      Cp_approx=Ggam_min(j)*Ug_fun(k_half,D_min,omAdd)+...
%                 Ggam_plus(j)*Ug_fun(k_half,D_plus,omAdd)+...
%                 Ggam_mid(j)*Ug_fun(k_half,D_mid,omAdd); 
%      err(j,:)=abs(Cp_exact-Cp_approx)./abs(Cp_exact);
%      indknu  =funC_closest(k_half,kappanu_D(j));
%      errknu(j)=abs(Cp_exact(indknu)-Cp_approx(indknu))/abs(Cp_exact(indknu));
%      Cnuv(j)=sqrt(Cp_exact(indknu));
%  end
%  
% subplot(4,1,2)
% mesh(DDepth,k_half,err');
% view(2);
% xlabel('depth (m) ');
% ylabel('wave number');
% axis([min(DDepth) max(DDepth) 0 max(k_half)]);
% hold on;   
% h=colorbar;
% %caxis([0 1]);
% %title('Error of Interpolation of Cp')
% ylabel(h, '|Exact-Approx|/|Exact|');
% lineobj = findobj('type', 'line');
% set(lineobj, 'linewidth', 1.3);
% textobj = findobj('type', 'text');
% set(textobj, 'fontunits', 'points');
% set(textobj, 'fontsize', 12);
% set(textobj, 'fontweight', 'bold');
% %saveas(gcf,['ErrIn3p'],'fig');
% plot_properties;
% set(textobj, 'fontsize', 18);
% set(lineobj, 'linewidth', 2.3);
%  %Error at wavenumber in the peak frequency
% subplot(4,1,3)
% plot(DDepth,errknu,'r');
% %title('Error of Cp (red) at peak frequencies')
% xlabel('depth (m)');
% ylabel('|Exact-Approx|/|Exact|');
% xlim([min(DDepth) max(DDepth)]);
% grid on;    
% lineobj = findobj('type', 'line');
% set(lineobj, 'linewidth', 1.3);
% textobj = findobj('type', 'text');
% set(textobj, 'fontunits', 'points');
% set(textobj, 'fontsize', 12);
% set(textobj, 'fontweight', 'bold');
% plot_properties;
% 
% set(textobj, 'fontsize', 18);
% set(lineobj, 'linewidth', 2.3);
% %saveas(gcf,['errknu'],'fig');          
% 
% subplot(4,1,4)
% plot(DDepth,kappanu_D);
% title('wave number at peak frequency');
% xlabel('depth');
% ylabel('\kappa');