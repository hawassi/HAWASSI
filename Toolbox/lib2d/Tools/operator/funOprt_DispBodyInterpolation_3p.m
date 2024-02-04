function [InterpDisp]= funOprt_DispBodyInterpolation_3p(g,dom,Ug_fun,Up_fun,Om_fun,omAdd,nupeak,Depth,bath,...
    Proj,dynmodel)
% 2-parameter approx of dispersion relation by interpolation
% cond: Cp  and Om^2 exact at kapnu for d Om^2 operator

nu      = mean(nupeak);
if strcmp(bath.type,'Sx');
    k_cut=max(dom.kx)/dom.cutfracx ;
    D_min=(nu/k_cut)^2/g;
elseif strcmp(bath.type,'Sy');
    k_cut=max(dom.ky)/dom.cutfracy ;
    D_min=(nu/k_cut)^2/g;
else
    D_min   = min(min(Depth));
end
D_plus  = max(max(Depth));
D_mid   = (D_plus+D_min)/2;
DDepth   =linspace(D_min,D_plus,100);

kappanu_D   = zeros(length(DDepth),1);
OneperOmSq_D        = zeros(length(DDepth),1);
OneperOmSq_min      = zeros(length(DDepth),1);
OneperOmSq_mid      = zeros(length(DDepth),1);
OneperOmSq_plus     = zeros(length(DDepth),1);
OmSq_D        = zeros(length(DDepth),1);
OmSq_min      = zeros(length(DDepth),1);
OmSq_mid      = zeros(length(DDepth),1);
OmSq_plus     = zeros(length(DDepth),1);
Cp_D        = zeros(length(DDepth),1);
Cp_min      = zeros(length(DDepth),1);
Cp_mid      = zeros(length(DDepth),1);
Cp_plus     = zeros(length(DDepth),1);

% if strcmp(bath.type,'Bx')
%     k=dom.kx;
% else
%     k=dom.ky;
% end
kk=sort(reshape(sqrt(dom.Kx.^2+dom.Ky.^2),[],1),'ascend');
k=linspace(min(kk),max(kk),200);

for j= 1:length(DDepth)
    D               = DDepth(j);
    kappanu_D(j)    = interp1(Om_fun(k,D,g),k,nu);%interpolation to find k at peak frequency nu
    OneperOmSq_min(j) = 1./Om_fun(kappanu_D(j),D_min,g,omAdd).^2;%sqrt(9.81*D_min);%Ug_fun(kappanu_D(j),D_min,omAdd);
    OneperOmSq_mid(j) = 1./Om_fun(kappanu_D(j),D_mid,g,omAdd).^2;%sqrt(9.81*D_mid);%Ug_fun(kappanu_D(j),D_mid,omAdd);
    OneperOmSq_plus(j)= 1./Om_fun(kappanu_D(j),D_plus,g,omAdd).^2;%sqrt(9.81*D_plus);%Ug_fun(kappanu_D(j),D_plus,omAdd);
    OneperOmSq_D(j)   = 1./Om_fun(kappanu_D(j),D,g,omAdd).^2;%sqrt(9.81*D);%Ug_fun(kappanu_D(j),D,omAdd);
    OmSq_min(j)       = Om_fun(kappanu_D(j),D_min,g,omAdd).^2;
    OmSq_mid(j)       = Om_fun(kappanu_D(j),D_mid,g,omAdd).^2;
    OmSq_plus(j)      = Om_fun(kappanu_D(j),D_plus,g,omAdd).^2;
    OmSq_D(j)         = Om_fun(kappanu_D(j),D,g,omAdd).^2;
    Cp_min(j)       = Up_fun(kappanu_D(j),D_min,g,omAdd);
    Cp_mid(j)       = Up_fun(kappanu_D(j),D_mid,g,omAdd);
    Cp_plus(j)      = Up_fun(kappanu_D(j),D_plus,g,omAdd);
    Cp_D(j)         = Up_fun(kappanu_D(j),D,g,omAdd);
    
end

A= (OmSq_mid.*Cp_plus)-(OmSq_plus.*Cp_mid);
B=-(OmSq_min.*Cp_plus)+(OmSq_plus.*Cp_min);
C= (OmSq_min.*Cp_mid)-(OmSq_mid.*Cp_min);
D=-(OneperOmSq_mid.*Cp_plus)+(OneperOmSq_plus.*Cp_mid);
E= (OneperOmSq_min.*Cp_plus)-(OneperOmSq_plus.*Cp_min);

F=-(OneperOmSq_min.*Cp_mid)+(OneperOmSq_mid.*Cp_min);
G= (OneperOmSq_mid.*OmSq_plus)-(OneperOmSq_plus.*OmSq_mid);
H=-(OneperOmSq_min.*OmSq_plus)+(OneperOmSq_plus.*OmSq_min);
I= (OneperOmSq_min.*OmSq_mid)-(OneperOmSq_mid.*OmSq_min);


det=OneperOmSq_min.*A+OneperOmSq_mid.*B+OneperOmSq_plus.*C;


Ggam_min  = ( A.*OneperOmSq_D+D.*OmSq_D+G.*Cp_D )./det;
Ggam_mid  = ( B.*OneperOmSq_D+E.*OmSq_D+H.*Cp_D )./det;
Ggam_plus = ( C.*OneperOmSq_D+F.*OmSq_D+I.*Cp_D )./det;

gam_plus=interp1(DDepth,Ggam_plus,Depth,'spline');
gam_min=interp1(DDepth,Ggam_min,Depth,'spline');
gam_mid=interp1(DDepth,Ggam_mid,Depth,'spline');
InterpDisp.gam_plus=gam_plus;
InterpDisp.gam_min= gam_min;
InterpDisp.gam_mid= gam_mid;
InterpDisp.D_min  =D_min;
InterpDisp.D_mid  =D_mid;
InterpDisp.D_plus =D_plus;


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
figure
subplot(4,1,1)
plot(DDepth,Ggam_min,'r',...
    DDepth,Ggam_mid,'g',DDepth,Ggam_plus,'-.b');
title([Proj.savename,', ',dynmodel,', Interpol functions for minimal(b) and maximal depth vs depth']);
%saveas(gcf,[folder_root,savename,dynmodel,'gamminmidplus.fig'])

%  pause
k_half=k(1:floor(end/8))';
err   =zeros(length(DDepth),length(k_half));
errknu=zeros(length(DDepth),1);
Omnuv=zeros(length(DDepth),1);
for j=1:length(DDepth)
    Om2_exact=Om_fun(k_half,DDepth(j),g,omAdd).^2;
    Om2_approx=Ggam_min(j)*Om_fun(k_half,D_min,g,omAdd).^2+...
        Ggam_plus(j)*Om_fun(k_half,D_plus,g,omAdd).^2+...
        Ggam_mid(j)*Om_fun(k_half,D_mid,g,omAdd).^2;
    err(j,:)=abs(Om2_exact-Om2_approx)./abs(Om2_exact);
    indknu  =funC_closest(k_half,kappanu_D(j));
    errknu(j)=abs(Om2_exact(indknu)-Om2_approx(indknu))/abs(Om2_exact(indknu));
    Omnuv(j)=sqrt(Om2_exact(indknu));
end

subplot(4,1,2)
mesh(DDepth,k_half,err');
view(2);
xlabel('depth (m) ');
ylabel('wave number');
axis([min(DDepth) max(DDepth) 0 max(k_half)]);
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