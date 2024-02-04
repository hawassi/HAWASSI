function [sg,H_min1,H_plus1,H_mid1,H_min2,H_plus2,H_mid2]= Up3IP_2part_runup(k,Up,Om,omAdd,nupeak,cutfrac,g,Depth,H_min,H_mid,H_plus,...
    folder_root,savename,model,infcategory)
% 3-parameter approx of phase velocity by interpolation
% cond: C, C0 and dC/deta exact exact at kapnu for disp of  C, C* and dC/deta operator

% H_min=5;
Htot1 = linspace(H_min,H_mid,100);
Htot2 = linspace(H_mid,H_plus,100);
H_min1=H_min;
H_plus1=H_mid;
H_mid1=(H_min1+H_plus1)/2;
H_min2=H_mid;
H_plus2=H_plus;

if strcmp(infcategory,'Deep water')
    Ddeep = invUpCgDeep(nupeak,H_min2,H_plus/300,0.51);
    if Ddeep==H_min2
    Ddeep=H_plus2;
    end
    H_mid2=H_min2+(Ddeep-H_min2)/6; %H_min2+(H_plus2-H_min2)/12;
   
else
    H_mid2=H_min2+(H_plus2-H_min2)/2;
end
k_cut=k(floor(length(k)/cutfrac));
wcut=Om(k_cut,H_plus,omAdd);

if strcmp(model.dispersion,'OmSWE')
    sg=[];
else
    Ninterp=length(Htot1);
    Cp2_H1 = zeros(Ninterp,1);
    Cp2_min1 = zeros(Ninterp,1);
    Cp2_plus1 = zeros(Ninterp,1);
    Cp2_mid1 = zeros(Ninterp,1);
    Cp2kL_H1 = zeros(Ninterp,1);
    Cp2kL_min1 = zeros(Ninterp,1);
    Cp2kL_plus1 = zeros(Ninterp,1);
    Cp2kL_mid1 = zeros(Ninterp,1);
    Cp2IG_H1 = zeros(Ninterp,1);
    Cp2IG_min1 = zeros(Ninterp,1);
    Cp2IG_plus1 = zeros(Ninterp,1);
    Cp2IG_mid1 = zeros(Ninterp,1);
    
    Cp2_H2 = zeros(Ninterp,1);
    Cp2_min2 = zeros(Ninterp,1);
    Cp2_plus2 = zeros(Ninterp,1);
    Cp2_mid2 = zeros(Ninterp,1);
    Cp2kL_H2 = zeros(Ninterp,1);
    Cp2kL_min2 = zeros(Ninterp,1);
    Cp2kL_plus2 = zeros(Ninterp,1);
    Cp2kL_mid2 = zeros(Ninterp,1);
    Cp2IG_H2 = zeros(Ninterp,1);
    Cp2IG_min2 = zeros(Ninterp,1);
    Cp2IG_plus2 = zeros(Ninterp,1);
    Cp2IG_mid2 = zeros(Ninterp,1);
    
    k_half      =k(1:floor(end/2));
    knu1         = zeros(Ninterp,1);
    knu2         = zeros(Ninterp,1);
    
    for j= 1:Ninterp
        Hj1=Htot1(j);
        Hj2=Htot2(j);
        %kappa0_H= interp1(Om(k,Hj),k,k(1));
        if strcmp(model.dispersion,'OmKdV')||strcmp(model.dispersion,'OmBBM')
%             Omega       =Om(k_half,Hj1,omAdd);
%             [~,Idnu]    =max(Omega);
            kappanu_H1   = funOprt_invOmExactVal(nupeak,Hj1,g);%Rinterp1(Omega(1:Idnu),k_half(1:Idnu),nupeak,'spline');% invOmExact(2*pi/nupeak,depth);
        else
            
            
            kappanu_H1   = funOprt_invOmExactVal(nupeak,Hj1,g);%interp1(Om(k_half,Hj1,omAdd),k_half,nupeak,'spline');% invOmExact(2*pi/nupeak,depth);
            kappanuIG_H1 = funOprt_invOmExactVal(nupeak/10,Hj1,g);%interp1(Om(k_half,Hj1,omAdd),k_half,nupeak/10,'spline');
            kappanuL_H1  = funOprt_invOmExactVal(wcut,Hj1,g);%interp1(Om(k_half,Hj1,omAdd),k_half,wcut,'spline');
            kappanu_H2   = funOprt_invOmExactVal(nupeak,Hj2,g);% invOmExact(2*pi/nupeak,depth);
            kappanuIG_H2 = funOprt_invOmExactVal(nupeak/10,Hj2,g);%interp1(Om(k_half,Hj2,omAdd),k_half,nupeak/10,'spline');
            kappanuL_H2  = funOprt_invOmExactVal(wcut,Hj2,g);%interp1(Om(k_half,Hj2,omAdd),k_half,wcut,'spline');
        end
        knu1(j)=  kappanu_H1;knu2(j)=  kappanu_H2;
        %   if   nupeak<0.5
        %       kappanu_H   = interp1(Om(k_half,Hj,omAdd),k_half,20*nupeak,'spline');%
        %   end
        %kappa_n_nu_H= k_max;%interp1(Om(k,Hj),k,nu);
        %interpolation to find k at peak frequency nu
        
        Cp2_min1(j) = Up(kappanu_H1,H_min1,omAdd).^2;
        Cp2_plus1(j)= Up(kappanu_H1,H_plus1,omAdd).^2;
        Cp2_mid1(j) = Up(kappanu_H1,H_mid1,omAdd).^2;
        Cp2_H1(j)   = Up(kappanu_H1,Hj1,omAdd).^2;
        
        Cp2kL_min1(j) = Up(kappanuL_H1,H_min1,omAdd).^2;% sqrt(g*H_min) ;%
        Cp2kL_plus1(j)= Up(kappanuL_H1,H_plus1,omAdd).^2;%%sqrt(g*H_plus) ;%
        Cp2kL_mid1(j) = Up(kappanuL_H1,H_mid1,omAdd).^2;%sqrt(g*H_mid) ;%
        Cp2kL_H1(j)   = Up(kappanuL_H1,Hj1,omAdd).^2;% sqrt(g*Hj) ;%
        
        Cp2IG_min1(j) =  Up(kappanuIG_H1,H_min1,omAdd).^2;%
        Cp2IG_plus1(j)= Up(kappanuIG_H1,H_plus1,omAdd).^2;%%
        Cp2IG_mid1(j) = Up(kappanuIG_H1,H_mid1,omAdd).^2;%
        Cp2IG_H1(j)   = Up(kappanuIG_H1,Hj1,omAdd).^2;%
        
        Cp2_min2(j) = Up(kappanu_H2,H_min2,omAdd).^2;
        Cp2_plus2(j)= Up(kappanu_H2,H_plus2,omAdd).^2;
        Cp2_mid2(j) = Up(kappanu_H2,H_mid2,omAdd).^2;
        Cp2_H2(j)   = Up(kappanu_H2,Hj2,omAdd).^2;
        
        Cp2kL_min2(j) = Up(kappanuL_H2,H_min2,omAdd).^2;% sqrt(g*H_min) ;%
        Cp2kL_plus2(j)= Up(kappanuL_H2,H_plus2,omAdd).^2;%%sqrt(g*H_plus) ;%
        Cp2kL_mid2(j) = Up(kappanuL_H2,H_mid2,omAdd).^2;%sqrt(g*H_mid) ;%
        Cp2kL_H2(j)   = Up(kappanuL_H2,Hj2,omAdd).^2;% sqrt(g*Hj) ;%
        
        Cp2IG_min2(j) =  Up(kappanuIG_H2,H_min2,omAdd).^2;%
        Cp2IG_plus2(j)= Up(kappanuIG_H2,H_plus2,omAdd).^2;%%
        Cp2IG_mid2(j) = Up(kappanuIG_H2,H_mid2,omAdd).^2;%
        Cp2IG_H2(j)   = Up(kappanuIG_H2,Hj2,omAdd).^2;%
    end
    
    
    A1= (Cp2kL_mid1.*Cp2IG_plus1)-(Cp2kL_plus1.*Cp2IG_mid1);
    B1=-(Cp2kL_min1.*Cp2IG_plus1)+(Cp2kL_plus1.*Cp2IG_min1);
    C1= (Cp2kL_min1.*Cp2IG_mid1)-(Cp2kL_mid1.*Cp2IG_min1);
    D1=-(Cp2_mid1.*Cp2IG_plus1)+(Cp2_plus1.*Cp2IG_mid1);
    E1= (Cp2_min1.*Cp2IG_plus1)-(Cp2_plus1.*Cp2IG_min1);
    
    F1=-(Cp2_min1.*Cp2IG_mid1)+(Cp2_mid1.*Cp2IG_min1);
    G1= (Cp2_mid1.*Cp2kL_plus1)-(Cp2_plus1.*Cp2kL_mid1);
    H1=-(Cp2_min1.*Cp2kL_plus1)+(Cp2_plus1.*Cp2kL_min1);
    I1= (Cp2_min1.*Cp2kL_mid1)-(Cp2_mid1.*Cp2kL_min1);
    
    
    det1=Cp2_min1.*A1+Cp2_mid1.*B1+Cp2_plus1.*C1;
    
    gam_min1  = ( A1.*Cp2_H1+D1.*Cp2kL_H1+G1.*Cp2IG_H1 )./det1;
    gam_mid1  = ( B1.*Cp2_H1+E1.*Cp2kL_H1+H1.*Cp2IG_H1 )./det1;
    gam_plus1 = ( C1.*Cp2_H1+F1.*Cp2kL_H1+I1.*Cp2IG_H1 )./det1;
    
    A2= (Cp2kL_mid2.*Cp2IG_plus2)-(Cp2kL_plus2.*Cp2IG_mid2);
    B2=-(Cp2kL_min2.*Cp2IG_plus2)+(Cp2kL_plus2.*Cp2IG_min2);
    C2= (Cp2kL_min2.*Cp2IG_mid2)-(Cp2kL_mid2.*Cp2IG_min2);
    D2=-(Cp2_mid2.*Cp2IG_plus2)+(Cp2_plus2.*Cp2IG_mid2);
    E2= (Cp2_min2.*Cp2IG_plus2)-(Cp2_plus2.*Cp2IG_min2);
    
    F2=-(Cp2_min2.*Cp2IG_mid2)+(Cp2_mid2.*Cp2IG_min2);
    G2= (Cp2_mid2.*Cp2kL_plus2)-(Cp2_plus2.*Cp2kL_mid2);
    H2=-(Cp2_min2.*Cp2kL_plus2)+(Cp2_plus2.*Cp2kL_min2);
    I2= (Cp2_min2.*Cp2kL_mid2)-(Cp2_mid2.*Cp2kL_min2);
    
    
    det2=Cp2_min2.*A2+Cp2_mid2.*B2+Cp2_plus2.*C2;
    
    gam_min2  = ( A2.*Cp2_H2+D2.*Cp2kL_H2+G2.*Cp2IG_H2 )./det2;
    gam_mid2  = ( B2.*Cp2_H2+E2.*Cp2kL_H2+H2.*Cp2IG_H2 )./det2;
    gam_plus2 = ( C2.*Cp2_H2+F2.*Cp2kL_H2+I2.*Cp2IG_H2 )./det2;
    
    
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
%     if model.nonlinear==2   %% 200917 now second order uses also polynomial model
%         sg.p1=spline(HtotComb,gam_plus1tot);
%         sg.c1=spline(HtotComb,gam_mid1tot);
%         sg.m1=spline(HtotComb,gam_min1tot);
%         sg.p2=spline(HtotComb,gam_plus2tot);
%         sg.c2=spline(HtotComb,gam_mid2tot);
%         sg.m2=spline(HtotComb,gam_min2tot);
%         %  gam_plus1tot=ppval(sg.p1,HtotComb);
%         %  gam_mid1tot=ppval(sg.c1,HtotComb);
%         %  gam_min1tot=ppval(sg.m1,HtotComb);
%         %  gam_plus2tot=ppval(sg.p2,HtotComb);
%         %  gam_mid2tot=ppval(sg.c2,HtotComb);
%         %  gam_min2tot=ppval(sg.m2,HtotComb);
%         
%     else

        DepthI=Depth;
        sg.p1=interp1(HtotComb,gam_plus1tot,DepthI,'spline');
        sg.c1=interp1(HtotComb,gam_mid1tot,DepthI,'spline');
        sg.m1=interp1(HtotComb,gam_min1tot,DepthI,'spline');
        sg.p2=interp1(HtotComb,gam_plus2tot,DepthI,'spline');
        sg.c2=interp1(HtotComb,gam_mid2tot,DepthI,'spline');
        sg.m2=interp1(HtotComb,gam_min2tot,DepthI,'spline');
        sg.p1(Depth<H_min)=0;sg.c1(Depth<H_min)=0;sg.m1(Depth<H_min)=0;
        sg.p2(Depth<H_min)=0;sg.c2(Depth<H_min)=0;sg.m2(Depth<H_min)=0;
        
        
        sg.Sp_p1=spline(HtotComb,gam_plus1tot);
        sg.Sp_c1=spline(HtotComb,gam_mid1tot);
        sg.Sp_m1=spline(HtotComb,gam_min1tot);
        sg.Sp_p2=spline(HtotComb,gam_plus2tot);
        sg.Sp_c2=spline(HtotComb,gam_mid2tot);
        sg.Sp_m2=spline(HtotComb,gam_min2tot);
       
%          figure;
%          plot(Depth,sg.p1,'r',Depth,sg.c1,'g',Depth,sg.m1,'b')
%    end
    
    
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
    
    knuj=zeros(size(Htotnow));
    for j=1:length(Htotnow)
        Cp_exact=Up(k_half,Htotnow(j),omAdd).^2;
        if j<indHmin
            Cp_approx=(sqrt(g*Htotnow(j)).*ones(size(k_half))).^2;
        elseif j<indHmid
            Cp_approx=gam_min1(j-indHmin+1)*Up(k_half,H_min1,omAdd).^2+...
                gam_plus1(j-indHmin+1)*Up(k_half,H_plus1,omAdd).^2+...
                gam_mid1(j-indHmin+1)*Up(k_half,H_mid1,omAdd).^2;
        else
            Cp_approx=gam_min2(j-indHmid+1).*Up(k_half,H_min2,omAdd).^2+...
                gam_plus2(j-indHmid+1).*Up(k_half,H_plus2,omAdd).^2+...
                gam_mid2(j-indHmid+1).*Up(k_half,H_mid2,omAdd).^2;
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
    axis([0 max(Htotnow) 0 k_cut]);
    hold on;
    h=colorbar;
    indKcut=closest(k_half,k_cut);
    
    caxis([0 max(max(err(:,1:indKcut)))]);
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
