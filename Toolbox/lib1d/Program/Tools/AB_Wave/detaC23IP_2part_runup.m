function [sg]= detaC23IP_2part_runup(k,Up,Om,omAdd,detaC2,nupeak,g,Depth,H_min,H_mid,...
        H_plus,folder_root,savename,model)
% 3-parameter approx of phase velocity by interpolation
% cond: C, C0 and dC/deta exact exact at kapnu for disp of  C, C* and dC/deta operator
    

Htot1 = linspace(H_min,H_mid,100);
Htot2 = linspace(H_mid,H_plus,100);
H_min1=H_min;
H_plus1=H_mid;
H_mid1=(H_min1+H_plus1)/2;
H_min2=H_mid;
H_plus2=H_plus;
H_mid2=(H_min2+H_plus2)/2;



if strcmp(model.dispersion,'OmSWE')
sg=[];
else
Ninterp=length(Htot1);
detadetaC22_H1 = zeros(Ninterp,1);
detadetaC22_min1 = zeros(Ninterp,1);
detadetaC22_plus1 = zeros(Ninterp,1);
detadetaC22_mid1 = zeros(Ninterp,1);
detadetaC22kL_H1 = zeros(Ninterp,1);
detadetaC22kL_min1 = zeros(Ninterp,1);
detadetaC22kL_plus1 = zeros(Ninterp,1);
detadetaC22kL_mid1 = zeros(Ninterp,1);
detadetaC22IG_H1 = zeros(Ninterp,1);
detadetaC22IG_min1 = zeros(Ninterp,1);
detadetaC22IG_plus1 = zeros(Ninterp,1);
detadetaC22IG_mid1 = zeros(Ninterp,1);

detadetaC22_H2 = zeros(Ninterp,1);
detadetaC22_min2 = zeros(Ninterp,1);
detadetaC22_plus2 = zeros(Ninterp,1);
detadetaC22_mid2 = zeros(Ninterp,1);
detadetaC22kL_H2 = zeros(Ninterp,1);
detadetaC22kL_min2 = zeros(Ninterp,1);
detadetaC22kL_plus2 = zeros(Ninterp,1);
detadetaC22kL_mid2 = zeros(Ninterp,1);
detadetaC22IG_H2 = zeros(Ninterp,1);
detadetaC22IG_min2 = zeros(Ninterp,1);
detadetaC22IG_plus2 = zeros(Ninterp,1);
detadetaC22IG_mid2 = zeros(Ninterp,1);

k_half      =k(1:floor(end/2));
knu1         = zeros(Ninterp,1);
knu2         = zeros(Ninterp,1);

for j= 1:Ninterp
Hj1=Htot1(j);
Hj2=Htot2(j);
%kappa0_H= interp1(Om(k,Hj),k,k(1));
if strcmp(model.dispersion,'OmKdV')||strcmp(model.dispersion,'OmBBM')
    Omega       =Om(k_half,Hj1,omAdd);
    [~,Idnu]    =max(Omega);
    kappanu_H1   = interp1(Omega(1:Idnu),k_half(1:Idnu),nupeak,'spline');% invOmExact(2*pi/nupeak,depth);
else
    kappanu_H1   = interp1(Om(k_half,Hj1,omAdd),k_half,nupeak,'spline');% invOmExact(2*pi/nupeak,depth);
    kappanuIG_H1   = interp1(Om(k_half,Hj1,omAdd),k_half,nupeak/10,'spline');
    kappanuL_H1 = interp1(Om(k_half,Hj1,omAdd),k_half,nupeak*6,'spline');
    kappanu_H2   = interp1(Om(k_half,Hj2,omAdd),k_half,nupeak,'spline');% invOmExact(2*pi/nupeak,depth);
    kappanuIG_H2   = interp1(Om(k_half,Hj2,omAdd),k_half,nupeak/10,'spline');
    kappanuL_H2 = interp1(Om(k_half,Hj2,omAdd),k_half,nupeak*3,'spline');
end
 knu1(j)=  kappanu_H1;knu2(j)=  kappanu_H2;
%   if   nupeak<0.5
%       kappanu_H   = interp1(Om(k_half,Hj,omAdd),k_half,20*nupeak,'spline');%
%   end
%kappa_n_nu_H= k_max;%interp1(Om(k,Hj),k,nu);
%interpolation to find k at peak frequency nu

detadetaC22_min1(j) = detaC2(kappanu_H1,H_min1);
detadetaC22_plus1(j)= detaC2(kappanu_H1,H_plus1);
detadetaC22_mid1(j) = detaC2(kappanu_H1,H_mid1);
detadetaC22_H1(j)   = detaC2(kappanu_H1,Hj1);

detadetaC22kL_min1(j) = detaC2(kappanuL_H1,H_min1);% sqrt(g*H_min) ;%
detadetaC22kL_plus1(j)= detaC2(kappanuL_H1,H_plus1);%%sqrt(g*H_plus) ;%
detadetaC22kL_mid1(j) = detaC2(kappanuL_H1,H_mid1);%sqrt(g*H_mid) ;%
detadetaC22kL_H1(j)   = detaC2(kappanuL_H1,Hj1);% sqrt(g*Hj) ;%

detadetaC22IG_min1(j) =  detaC2(kappanuIG_H1,H_min1);%
detadetaC22IG_plus1(j)= detaC2(kappanuIG_H1,H_plus1);%%
detadetaC22IG_mid1(j) = detaC2(kappanuIG_H1,H_mid1);%
detadetaC22IG_H1(j)   = detaC2(kappanuIG_H1,Hj1);%

detadetaC22_min2(j) = detaC2(kappanu_H2,H_min2);
detadetaC22_plus2(j)= detaC2(kappanu_H2,H_plus2);
detadetaC22_mid2(j) = detaC2(kappanu_H2,H_mid2);
detadetaC22_H2(j)   = detaC2(kappanu_H2,Hj2);

detadetaC22kL_min2(j) = detaC2(kappanuL_H2,H_min2);% sqrt(g*H_min) ;%
detadetaC22kL_plus2(j)= detaC2(kappanuL_H2,H_plus2);%%sqrt(g*H_plus) ;%
detadetaC22kL_mid2(j) = detaC2(kappanuL_H2,H_mid2);%sqrt(g*H_mid) ;%
detadetaC22kL_H2(j)   = detaC2(kappanuL_H2,Hj2);% sqrt(g*Hj) ;%

detadetaC22IG_min2(j) =  detaC2(kappanuIG_H2,H_min2);%
detadetaC22IG_plus2(j)= detaC2(kappanuIG_H2,H_plus2);%%
detadetaC22IG_mid2(j) = detaC2(kappanuIG_H2,H_mid2);%
detadetaC22IG_H2(j)   = detaC2(kappanuIG_H2,Hj2);%
end


A1= (detadetaC22kL_mid1.*detadetaC22IG_plus1)-(detadetaC22kL_plus1.*detadetaC22IG_mid1);
B1=-(detadetaC22kL_min1.*detadetaC22IG_plus1)+(detadetaC22kL_plus1.*detadetaC22IG_min1);
C1= (detadetaC22kL_min1.*detadetaC22IG_mid1)-(detadetaC22kL_mid1.*detadetaC22IG_min1);
D1=-(detadetaC22_mid1.*detadetaC22IG_plus1)+(detadetaC22_plus1.*detadetaC22IG_mid1);
E1= (detadetaC22_min1.*detadetaC22IG_plus1)-(detadetaC22_plus1.*detadetaC22IG_min1);

F1=-(detadetaC22_min1.*detadetaC22IG_mid1)+(detadetaC22_mid1.*detadetaC22IG_min1);
G1= (detadetaC22_mid1.*detadetaC22kL_plus1)-(detadetaC22_plus1.*detadetaC22kL_mid1);
H1=-(detadetaC22_min1.*detadetaC22kL_plus1)+(detadetaC22_plus1.*detadetaC22kL_min1);
I1= (detadetaC22_min1.*detadetaC22kL_mid1)-(detadetaC22_mid1.*detadetaC22kL_min1);


det1=detadetaC22_min1.*A1+detadetaC22_mid1.*B1+detadetaC22_plus1.*C1;

gam_min1  = ( A1.*detadetaC22_H1+D1.*detadetaC22kL_H1+G1.*detadetaC22IG_H1 )./det1;
gam_mid1  = ( B1.*detadetaC22_H1+E1.*detadetaC22kL_H1+H1.*detadetaC22IG_H1 )./det1;
gam_plus1 = ( C1.*detadetaC22_H1+F1.*detadetaC22kL_H1+I1.*detadetaC22IG_H1 )./det1;

A2= (detadetaC22kL_mid2.*detadetaC22IG_plus2)-(detadetaC22kL_plus2.*detadetaC22IG_mid2);
B2=-(detadetaC22kL_min2.*detadetaC22IG_plus2)+(detadetaC22kL_plus2.*detadetaC22IG_min2);
C2= (detadetaC22kL_min2.*detadetaC22IG_mid2)-(detadetaC22kL_mid2.*detadetaC22IG_min2);
D2=-(detadetaC22_mid2.*detadetaC22IG_plus2)+(detadetaC22_plus2.*detadetaC22IG_mid2);
E2= (detadetaC22_min2.*detadetaC22IG_plus2)-(detadetaC22_plus2.*detadetaC22IG_min2);

F2=-(detadetaC22_min2.*detadetaC22IG_mid2)+(detadetaC22_mid2.*detadetaC22IG_min2);
G2= (detadetaC22_mid2.*detadetaC22kL_plus2)-(detadetaC22_plus2.*detadetaC22kL_mid2);
H2=-(detadetaC22_min2.*detadetaC22kL_plus2)+(detadetaC22_plus2.*detadetaC22kL_min2);
I2= (detadetaC22_min2.*detadetaC22kL_mid2)-(detadetaC22_mid2.*detadetaC22kL_min2);


det2=detadetaC22_min2.*A2+detadetaC22_mid2.*B2+detadetaC22_plus2.*C2;

gam_min2  = ( A2.*detadetaC22_H2+D2.*detadetaC22kL_H2+G2.*detadetaC22IG_H2 )./det2;
gam_mid2  = ( B2.*detadetaC22_H2+E2.*detadetaC22kL_H2+H2.*detadetaC22IG_H2 )./det2;
gam_plus2 = ( C2.*detadetaC22_H2+F2.*detadetaC22kL_H2+I2.*detadetaC22IG_H2 )./det2;


 if model.nonlinear==2
 HtotComb= [Htot1,Htot2(2:end)];  
 assignin('base','HtotComb',HtotComb)
 indHmin1=closest(HtotComb,H_min1);
 indHplus1=closest(HtotComb,H_plus1);
 indHplus2=closest(HtotComb,H_plus2);
 gam_min1tot=zeros(size(HtotComb));
 gam_mid1tot=zeros(size(HtotComb));
 gam_plus1tot=zeros(size(HtotComb));
 gam_min2tot=zeros(size(HtotComb));
 gam_mid2tot=zeros(size(HtotComb));
 gam_plus2tot=zeros(size(HtotComb));
 gam_min1tot(indHmin1:indHplus1)=gam_min1;
 gam_mid1tot(indHmin1:indHplus1)=gam_mid1;
 gam_plus1tot(indHmin1:indHplus1)=gam_plus1;
 gam_min2tot(indHplus1+1:indHplus2)=gam_min2(2:end).';
 gam_mid2tot(indHplus1+1:indHplus2)=gam_mid2(2:end).';
 gam_plus2tot(indHplus1+1:indHplus2)=gam_plus2(2:end).';
  
 sg.p1=spline(HtotComb,gam_plus1tot);
 sg.c1=spline(HtotComb,gam_mid1tot);
 sg.m1=spline(HtotComb,gam_min1tot);
 sg.p2=spline(HtotComb,gam_plus2tot);
 sg.c2=spline(HtotComb,gam_mid2tot);
 sg.m2=spline(HtotComb,gam_min2tot);
%  gam_plus1tot=ppval(sg.p1,HtotComb);
%  gam_mid1tot=ppval(sg.c1,HtotComb);
%  gam_min1tot=ppval(sg.m1,HtotComb);
%  gam_plus2tot=ppval(sg.p2,HtotComb);
%  gam_mid2tot=ppval(sg.c2,HtotComb);
%  gam_min2tot=ppval(sg.m2,HtotComb);
 
 else
  DepthI=Depth;DepthI(Depth<H_min)=H_min;  
 sg.p=interp1(Htot,gam_plus1,DepthI,'spline');
 sg.c=interp1(Htot,gam_mid1,DepthI,'spline');
 sg.m=interp1(Htot,gam_min1,DepthI,'spline');
%  figure;
%  plot(Depth,sg.p,'r',Depth,sg.c,'g',Depth,sg.m,'b')
 end

 hfigInterp=figure('Name','Bathy Set-up', 'Position',[170,10,2,600]);
 set(hfigInterp,'Renderer','zbuffer'); %due to graphics driver
 subplot(4,1,1)
%  plot([Htot1,Htot2],[gam_min1.',gam_min2.'],'r',...
%      [Htot1,Htot2],[gam_mid1.',gam_mid2.'],'--g',...
%      [Htot1,Htot2],[gam_plus1.',gam_plus2.'],'-.b');
  plot(HtotComb,gam_min1tot,'--r',HtotComb,gam_mid1tot,'--g',HtotComb,gam_plus1tot,'--b',...
      HtotComb,gam_min2tot,'r',HtotComb,gam_mid2tot,'g',HtotComb,gam_plus2tot,'b');
  title([savename,', ',model.dyn,', Interpol functions for minimal(b) and maximal depth vs depth']);
 %saveas(gcf,[folder_root,savename,dynmodel,'gamminmidplus.fig'])

 
%  pause
Htotnow=[linspace(0,H_min1-H_min1/10,20), Htot1, Htot2(2:end)];
assignin('base','Htotnow',Htotnow)
err   =zeros(length(Htotnow),length(k_half));
errknu=zeros(length(Htotnow),1);

indHmin=closest(Htotnow,H_min1);
indHmid=closest(Htotnow,H_plus1);

detaC2_exactN=zeros(length(k_half),length(Htotnow));
detaC2_approxN=zeros(length(k_half),length(Htotnow));

knuj=zeros(size(Htotnow));
for j=1:length(Htotnow)
     detaC2_exact=detaC2(k_half,Htotnow(j));
     if j<indHmin
     detaC2_approx=g.*ones(size(k_half));
     elseif j<indHmid
     detaC2_approx=gam_min1(j-indHmin+1)*detaC2(k_half,H_min1)+...
                gam_plus1(j-indHmin+1)*detaC2(k_half,H_plus1)+...
                gam_mid1(j-indHmin+1)*detaC2(k_half,H_mid1);     
     else
     detaC2_approx=gam_min2(j-indHmid+1).*detaC2(k_half,H_min2)+...
                gam_plus2(j-indHmid+1).*detaC2(k_half,H_plus2)+...
                gam_mid2(j-indHmid+1).*detaC2(k_half,H_mid2); 
     end
     err(j,:)=abs(detaC2_exact-detaC2_approx)./detaC2_exact;
    
     if j>1
     knuj(j)= interp1(Om(k_half,Htotnow(j),omAdd),k_half,nupeak,'spline');% 
     else
     knuj(j)=0.00001;    
     end
     indknu  =closest(k_half,knuj(j));
     errknu(j)=abs(detaC2_exact(indknu)-detaC2_approx(indknu))/detaC2_exact(indknu);
     detaC2_exactN(:,j)=detaC2_exact.';
     detaC2_approxN(:,j)=detaC2_approx.';
     
 end
 
subplot(4,1,2)
%mesh(Htotnow,k_half,err');
mesh(Htotnow,k_half,detaC2_exactN);
view(2);
xlabel('depth (m) ');
ylabel('wave number');
axis([0 max(Htotnow) 0 max(k_half)]);
hold on;
h=colorbar;
caxis([0 10]);
title('detaC2 exact.')
ylabel(h, 'detaC2 exact.');
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

% plot(Htotnow,errknu,'r');
%plot(Htotnow,errknuCder,'r');
% title('Error of detaC2 (red) at peak wavenumbers')
% xlabel('depth (m)');
% ylabel('|Exct-Approx|/Exct');
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

% detaC2_H_app=gam_min.*detaC2_min+gam_plus.*detaC2_plus+gam_mid.*detaC2_mid;
% C0_H_app=gam_min.*C0_min+gam_plus.*C0_plus+gam_mid.*C0_mid;
% Cder_H_app=gam_min.*Cder_min+gam_plus.*Cder_plus+gam_mid.*Cder_mid;

% figure;
% subplot(3,1,1)
% plot(Htot,detaC2_min,'r',Htot,detaC2_mid,'g',Htot,detaC2_plus,'b',Htot,detaC2_H,'k',Htot,detaC2_H_app,'--c')
% xlabel('depth')
% legend('detaC2_{min}','detaC2_{mid}','detaC2_{plus}','detaC2_{Exact}','detaC2_{Exact1}')
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
% plot(k_half,detaC2(k_half,H_min,omAdd),'r',k_half,detaC2(k_half,H_mid,omAdd),'g',k_half,detaC2(k_half,H_plus,omAdd),'b')
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
% detaC2j   =detaC2(k_half,totDepthObs(ii),omAdd);
% plot(wwj,detaC2j,colP{ii});
% hold on;
% end
% hold off;
% xlabel('\omega [rad/s]');
% ylabel('detaC2');
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
% detaC2j   = detaC2(k_half,totDepthObs(ii),omAdd);%
% plot(k_half,detaC2j,colP{ii});
% hold on;
% end
% hold off;
% xlabel('k');
% ylabel('detaC2 ');
% legend(['D=',num2str(totDepthObs(1))],['D=',num2str(totDepthObs(2))],...
%     ['D=',num2str(totDepthObs(3))],['D=',num2str(totDepthObs(4))],...
%     ['D=',num2str(totDepthObs(5))])
% plot_properties;
