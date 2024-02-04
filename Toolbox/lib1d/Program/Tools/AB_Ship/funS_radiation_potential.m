function [psi,Ma,dampC,Lpsiev,Fpsi,FpsiEv,dxi_Krad]=funS_radiation_potential(shippar,x,xShip,Sdraft,sZ,chiXcZ0,chiWL,nu,Smotion,D,Nship,etaL,etaR,dxetaL,dxetaR)

dx=x(2)-x(1);
ii=1;
    Idmotion=Smotion(ii);
    Xc=xShip(ii,2);%
    indl=closest(x,xShip(ii,1));indr=closest(x,xShip(ii,3));
    DL=D(indl);
    DR=D(indr);
    dxD=gradient(D,dx);
    dxDL=dxD(indl);dxDR=dxD(indr);
    
    HL=DL+etaL;
    HR=DR+etaR;
    modes=shippar.Evmode;
    b=shippar.form.Slength(1)/2;
    KappaL=shippar.kappa(1:modes+1);
    KappaR=shippar.kappa(modes+2:end);
    draft=Sdraft(ii)-sZ(ii);
   
    Z0=shippar.form.XZrot0(ii,2);%% for pitch
    X0=shippar.form.XZrot0(ii,1);%% always 0
    [Amat, BvectS, BvectH, BvectP] =...
        fun_construct_matrices_radiation(DL,DR,dxDL,dxDR,etaL,etaR,dxetaL,dxetaR,draft,b,...
        modes,KappaL,KappaR,Idmotion,Z0);
        
      
    if strcmpi(Idmotion,'Surge')|| strcmpi(Idmotion,'Free')
        CoefS=Amat\BvectS;
        ALS=CoefS(1);ARS=CoefS(2);tau0S=CoefS(3);tau1S=CoefS(4);
    end
    
    if strcmpi(Idmotion,'Heave')|| strcmpi(Idmotion,'Free')
        CoefH=Amat\BvectH;
        ALH=CoefH(1);ARH=CoefH(2);tau0H=CoefH(3);tau1H=CoefH(4);  
    end
    
    if strcmpi(Idmotion,'Pitch')|| strcmpi(Idmotion,'Free')
        CoefP=Amat\BvectP;
        ALP=CoefP(1);ARP=CoefP(2);tau0P=CoefP(3);tau1P=CoefP(4);
    end
     
%     assignin('base','Amat',Amat);
% %     assignin('base','BvectS',BvectS);
% %    assignin('base','BvectH',BvectH);
%      assignin('base','BvectP',BvectP);
% %     assignin('base','CoefS',CoefS);
% %    assignin('base','CoefH',CoefH);
%      assignin('base','CoefP',CoefP);
      
  
    %%% shifted to ship position xc;
    chiship=chiXcZ0(:,end)';
    chiSWL=chiship+chiWL(:,end)';
    indl=closest(x,Xc-b);indr=closest(x,Xc+b);
    
    FpsiEvS=zeros(size(x));FevLS=0;FevRS=0;
    FpsiS=zeros(size(x));
    FpsiEvH=zeros(size(x));FevLH=0;FevRH=0;
    FpsiH=zeros(size(x));
    FpsiEvP=zeros(size(x));FevLP=0;FevRP=0;
    FpsiP=zeros(size(x));
    
    %      Fpsi(indl)=AL*tan(kappa0*D);
    %      Fpsi(indr)=-BR*tan(kappa0*D);
    
    psiS=zeros(size(x));
    psiH=zeros(size(x));
    psiP=zeros(size(x));
    LpsievS=zeros(size(psiS));
    LpsievH=zeros(size(psiH));
    LpsievP=zeros(size(psiP));
    
    kappa0L=KappaL(1);
    kappaL=KappaL(2:end);
    kappa0R=KappaR(1);
    kappaR=KappaR(2:end);
    
    IntQ0TetaL=fun_IntegQn(kappa0L,etaL,DL,-draft,etaL);
    IntQ0TetaR=fun_IntegQn(kappa0R,etaR,DR,-draft,etaR);
    %%%for calculating d\xi K(\phi,\xi)
     ZshipL=linspace(-draft,etaL,100);
     ZshipR=linspace(-draft,etaR,100);
     DZL=DL+ZshipL;
     DZR=DR+ZshipR;
     dzQLz0=-kappa0L.*sin(kappa0L.*DZL)./cos(kappa0L.*HL);
     dzQRz0=-kappa0R.*sin(kappa0R.*DZR)./cos(kappa0R.*HR);
  
    
    if strcmpi(Idmotion,'Surge')|| strcmpi(Idmotion,'Free')
        f11=ALS.*IntQ0TetaL-ARS*IntQ0TetaR;
        psiS=real(tau1S.*(x-Xc)+tau0S).*chiSWL;
        
        dPhi1Sdz_zL=real(ALS.*dzQLz0);dPhi1Sdx_zL=1;
        dPhi3Sdz_zR=real(ARS.*dzQRz0);dPhi3Sdx_zR=1;
    end
    if strcmpi(Idmotion,'Heave')|| strcmpi(Idmotion,'Free')
        psiH=(((D-draft).^2-(x-Xc).^2)./2./(D-draft)+tau1H.*(x-Xc)+tau0H).*chiSWL;
        f33=trapz(x(indl:indr)-Xc,psiH(indl:indr));
        psiH=real(psiH);
        %%%for calculating d\xi K(\phi,\xi)
        dxGamma3_x=-(x-Xc)./(D-draft)+dxD.*(1-((D-draft)-(x-Xc).^2)./2./(D-draft));
        dxPsiH_x=real(dxGamma3_x+tau1H);
        dzPsiH_x=1 ;
        %%%%%%%%%%%
    end
    if strcmpi(Idmotion,'Pitch')|| strcmpi(Idmotion,'Free')
        psiP=(((D-draft).^2.*(x-Xc)-(x-Xc).^3/3)./2./(D-draft)+tau1P.*(x-Xc)+tau0P).*chiSWL;
         
       
        IZQ0TetaL = fun_IntegZQn(kappa0L,etaL,DL,-draft,etaL);
        IZQ0TetaR = fun_IntegZQn(kappa0R,etaR,DR,-draft,etaR);
        
        f551=ARP.*(IZQ0TetaR-Z0.*IntQ0TetaR)-ALP.*(IZQ0TetaL-Z0.*IntQ0TetaL);%
        f552=trapz(x(indl:indr)-Xc,(x(indl:indr)-Xc).*psiP(indl:indr));
         psiP=real(psiP);
        %%%for calculating d\xi K(\phi,\xi)
        dxGamma5_x=((D-draft).^2-(x-Xc).^2)./2./(D-draft)+dxD.*((x-Xc)-((D-draft).*(x-Xc)-(x-Xc).^3/3)./2./(D-draft));
        dxPsiP_x=real(dxGamma5_x+tau1P);
        dzPsiP_x=1 ;
        dPhi1Pdz_zL=real(ALP.*dzQLz0);dPhi1Pdx_zL=1;
        dPhi3Pdz_zR=real(ARP.*dzQRz0);dPhi3Pdx_zR=1;
        
    end
    
    if strcmpi(Idmotion,'Free')    
          f13= 2*tau0S*b;
          f151=ARS.*(IZQ0TetaR-Z0.*IntQ0TetaR)-ALS.*(IZQ0TetaL-Z0.*IntQ0TetaL);%
          f152=(-tau0S*(X0*2*b)+tau1S*2*b^3/3);
          f31=ALH.*IntQ0TetaL-ARH*IntQ0TetaR;   
          f351=ARH.*(IZQ0TetaR-Z0.*IntQ0TetaR)-ALH.*(IZQ0TetaL-Z0.*IntQ0TetaL);%;
          f352=trapz(x(indl:indr)-Xc,(x(indl:indr)-Xc).*psiH(indl:indr));
          f51=ALP.*IntQ0TetaL-ARP*IntQ0TetaR;
          f53= trapz(x(indl:indr)-Xc,psiP(indl:indr));
    end
    
    if modes>0
        for m=1:modes
            mu_m=m*pi./(D-draft);
            
             IntQmTetaL=fun_IntegQn(kappaL(m),etaL,DL,-draft,etaL);
             IntQmTetaR=fun_IntegQn(kappaR(m),etaR,DR,-draft,etaR);
                
            if strcmpi(Idmotion,'Surge')|| strcmpi(Idmotion,'Free')
                epsS1m=CoefS(4*m+1);epsS21m=CoefS(4*m+2);
                epsS22m=CoefS(4*m+3);epsS3m=CoefS(4*m+4);
                psiSev=real((epsS21m.*exp(mu_m.*((x-Xc)+b))+epsS22m.*exp(-mu_m.*((x-Xc)-b))).*cos(mu_m.*(D-draft))).*chiSWL;
                psiS=psiS+psiSev;
                
                f11=f11+(epsS1m.*IntQmTetaL-epsS3m.*IntQmTetaR);
                
                LpsievS(1:indl)=LpsievS(1:indl)-epsS1m*exp(kappaL(m).*(x(1:indl)-Xc+b)).*kappaL(m).*tan(kappaL(m).*D(1:indl));
                LpsievS(indr:end)=LpsievS(indr:end)-epsS3m*exp(-kappaR(m).*(x(indr:end)-Xc-b)).*kappaR(m).*tan(kappaR(m).*D(indr:end));
                FevLS=FevLS+epsS1m*tan(kappaL(m).*HL);
                FevRS=FevRS-epsS3m*tan(kappaR(m).*HR);
                
                %%%for calculating d\xi K(\phi,\xi)
                dzQLzm=-kappaL(m).*sin(kappaL(m).*DZL)./cos(kappaL(m).*HL);
                dzQRzm=-kappaR(m).*sin(kappaR(m).*DZR)./cos(kappaR(m).*HR);
   
                dPhi1Sdz_zL=dPhi1Sdz_zL+epsS1m.*dzQLzm;
                dPhi3Sdz_zR=dPhi3Sdz_zR+epsS3m.*dzQRzm;
            end
            if strcmpi(Idmotion,'Heave')|| strcmpi(Idmotion,'Free')
                epsH1m=CoefH(4*m+1);epsH21m=CoefH(4*m+2);
                epsH22m=CoefH(4*m+3);epsH3m=CoefH(4*m+4);
                psiHev=real((epsH21m.*exp(mu_m.*((x-Xc)+b))+epsH22m.*exp(-mu_m.*((x-Xc)-b))).*cos(mu_m.*(D-draft))).*chiSWL;
                psiH=psiH+psiHev;
                
              
                f33=f33+trapz(x(indl:indr)-Xc,psiHev(indl:indr));
                
                LpsievH(1:indl)=LpsievH(1:indl)-epsH1m*exp(kappaL(m).*(x(1:indl)-Xc+b)).*kappaL(m).*tan(kappaL(m).*D(1:indl));
                LpsievH(indr:end)=LpsievH(indr:end)-epsH3m*exp(-kappaR(m).*(x(indr:end)-Xc-b))*kappaR(m).*tan(kappaR(m).*D(indr:end));
                FevLH=FevLH+epsH1m*tan(kappaL(m)*HL);
                FevRH=FevRH-epsH3m*tan(kappaR(m)*HR);
                
                dxPsiHev=mu_m.*(1-(x-Xc).*dxD./(D-draft)).*(epsH21m.*exp(mu_m.*((x-Xc)+b))-epsH22m.*exp(-mu_m.*((x-Xc)-b))).*cos(mu_m.*(D-draft));
                dxPsiH_x=dxPsiH_x+real(dxPsiHev);
            end
            if strcmpi(Idmotion,'Pitch')|| strcmpi(Idmotion,'Free')
                epsP1m=CoefP(4*m+1);epsP21m=CoefP(4*m+2);
                epsP22m=CoefP(4*m+3);epsP3m=CoefP(4*m+4);
                psiPev=(epsP21m.*exp(mu_m.*((x-Xc)+b))+epsP22m.*exp(-mu_m.*((x-Xc)-b))).*cos(mu_m.*(D-draft)).*chiSWL;
                psiP=psiP+real(psiPev);
                
                IZQmTetaL = fun_IntegZQn(kappaL(m),etaL,DL,-draft,etaL);
                IZQmTetaR = fun_IntegZQn(kappaR(m),etaR,DR,-draft,etaR);
                
                f551=f551+epsP3m.*(IZQmTetaR-Z0.*IntQmTetaR)-epsP1m.*(IZQmTetaL-Z0.*IntQmTetaL);%
                f552=f552+trapz(x(indl:indr)-Xc,(x(indl:indr)-Xc).*psiPev(indl:indr));
                
                
                LpsievP(1:indl)=LpsievP(1:indl)-epsP1m*exp(kappaL(m).*(x(1:indl)-Xc+b)).*kappaL(m).*tan(kappaL(m).*D(1:indl));
                LpsievP(indr:end)=LpsievP(indr:end)-epsP3m.*exp(-kappaR(m).*(x(indr:end)-Xc-b)).*kappaR(m).*tan(kappaR(m).*D(indr:end));
                FevLP=FevLP+epsP1m*tan(kappaL(m)*HL);
                FevRP=FevRP-epsP3m*tan(kappaR(m)*HR);
                
                
                 %%%for calculating d\xi K(\phi,\xi)
                dzQLzm=-kappaL(m).*sin(kappaL(m).*DZL)./cos(kappaL(m).*HL);
                dzQRzm=-kappaR(m).*sin(kappaR(m).*DZR)./cos(kappaR(m).*HR);
   
                dPhi1Pdz_zL=dPhi1Pdz_zL+epsP1m.*dzQLzm;
                dPhi3Pdz_zR=dPhi3Pdz_zR+epsP3m.*dzQRzm;
                
                dxPsiPev=mu_m.*(1-(x-Xc).*dxD./(D-draft)).*(epsP21m.*exp(mu_m.*((x-Xc)+b))-epsP22m.*exp(-mu_m.*((x-Xc)-b))).*cos(mu_m.*(D-draft));
                dxPsiP_x=dxPsiP_x+real(dxPsiPev);
                
            end
            if strcmpi(Idmotion,'Free')
                 f13=f13+trapz(x(indl:indr)-Xc,psiSev(indl:indr));
                 f151=f151+epsS3m.*(IZQmTetaR-Z0.*IntQmTetaR)-epsS1m.*(IZQmTetaL-Z0.*IntQmTetaL);
                 f152=f152+trapz(x(indl:indr)-Xc,(x(indl:indr)-Xc).*psiSev(indl:indr));
                 f31=f31+(epsH1m.*IntQmTetaL-epsH3m.*IntQmTetaR);
                 f351=f351+epsH3m.*(IZQmTetaR-Z0.*IntQmTetaR)-epsH1m.*(IZQmTetaL-Z0.*IntQmTetaL);%
                 f352=f352+trapz(x(indl:indr)-Xc,(x(indl:indr)-Xc).*psiHev(indl:indr));
                 f51=f51+(epsP1m.*IntQmTetaL-epsP3m.*IntQmTetaR);
                 f53=f53+trapz(x(indl:indr)-Xc,psiPev(indl:indr));
            end
        end
    end
    %      figure;
    %      plot(x,psi);
    dxi_Krad.z(1)      = 0;
    dxi_Krad.x(1)      = 0;
    dxi_Krad.theta(1)   = 0;
    
    if strcmpi(Idmotion,'Surge')
        Ma33=0;Ma55=0;
        Ma11=real(f11);
        B11=imag(f11);B33=0;B55=0;
        dxi_Krad.x(ii)= real(trapz(ZshipL,(dPhi1Sdx_zL.^2+dPhi1Sdz_zL.^2)./2-dPhi1Sdz_zL.*(dPhi1Sdx_zL.*draft./dx)./sqrt(1+(draft./dx)^2)) ...
            -trapz(ZshipR,(dPhi3Sdx_zR.^2+dPhi3Sdz_zR.^2)./2+dPhi3Sdz_zR.*(dPhi3Sdx_zR.*draft./dx)./sqrt(1+(draft./dx)^2)));
    elseif strcmpi(Idmotion,'Heave')
        Ma33=real(f33);B33=imag(f33);
        Ma11=0;Ma55=0;B55=0;B11=0;
        % zMa=trapz(x(indl:indr),psi(indl:indr).*nutild.z(1,indl:indr));
        %      [zMa zMaN]
        
        dxi_Krad.z(1)     =  real(trapz(x(indl:indr)-Xc,(dxPsiH_x(indl:indr).^2+dzPsiH_x.^2)/2-dzPsiH_x.*nu.z(1,indl:indr)));%trapz(x(indl:indr)-Xc,(dxPsiH_x(indl:indr).^2-dzPsiH_x.^2)/2);%
    elseif strcmpi(Idmotion,'Pitch')
        f55=f551+f552;
        Ma55=real(f55);B55=imag(f55);B33=0;B11=0;
        %     thetaMaN=trapz(x(indl:indr),psi(indl:indr).*nutild.theta(1,indl:indr))
        Ma11=0;Ma33=0;
        dxi_Krad.theta(1)  = real(-trapz(ZshipL,(ZshipL-Z0).*((dPhi1Pdx_zL.^2+dPhi1Pdz_zL.^2)./2-dPhi1Pdz_zL.*(dPhi1Pdx_zL.*draft./dx)./sqrt(1+(draft./dx)^2))) ...
            +trapz(ZshipR,(ZshipL-Z0).*((dPhi3Pdx_zR.^2+dPhi3Pdz_zR.^2)./2+dPhi3Pdz_zR.*(dPhi3Pdx_zR.*draft./dx)./sqrt(1+(draft./dx)^2)))...
            +trapz(x(indl:indr)-Xc,(x(indl:indr)-X0-Xc).*((dxPsiP_x(indl:indr).^2+dzPsiP_x.^2)/2-dzPsiP_x.*nu.z(1,indl:indr))));      
    elseif strcmpi(Idmotion,'Free')
        f15=f151+f152;f35=f351+f352;f55=f551+f552;
        Ma11=real(f11);
        Ma13=real(f13);
        Ma15=real(f15);
        Ma31=real(f31);
        Ma33=real(f33); 
        Ma35=real(f35);
        Ma51=real(f51);
        Ma53=real(f53); 
        Ma55=real(f55);
      
         
        B11=imag(f11);B13=imag(f13); B15=imag(f15);
        B31=imag(f31);B33=imag(f33); B35=imag(f35);
        B51=imag(f51);B53=imag(f53); B55=imag(f55);
        dxi_Krad.z(1)     =  real(trapz(x(indl:indr)-Xc,(dxPsiH_x(indl:indr).^2+dzPsiH_x.^2)/2-dzPsiH_x.*nu.z(1,indl:indr)));%trapz(x(indl:indr)-Xc,(dxPsiH_x(indl:indr).^2-dzPsiH_x.^2)/2);%
        
        dxi_Krad.x(ii)= real(trapz(ZshipL,(dPhi1Sdx_zL.^2+dPhi1Sdz_zL.^2)./2-dPhi1Sdz_zL.*(dPhi1Sdx_zL.*draft./dx)./sqrt(1+(draft./dx)^2)) ...
            -trapz(ZshipR,(dPhi3Sdx_zR.^2+dPhi3Sdz_zR.^2)./2+dPhi3Sdz_zR.*(dPhi3Sdx_zR.*draft./dx)./sqrt(1+(draft./dx)^2)));
        
        dxi_Krad.theta(1)  = real(-trapz(ZshipL,(ZshipL-Z0).*((dPhi1Pdx_zL.^2+dPhi1Pdz_zL.^2)./2-dPhi1Pdz_zL.*(dPhi1Pdx_zL.*draft./dx)./sqrt(1+(draft./dx)^2))) ...
            +trapz(ZshipR,(ZshipL-Z0).*((dPhi3Pdx_zR.^2+dPhi3Pdz_zR.^2)./2+dPhi3Pdz_zR.*(dPhi3Pdx_zR.*draft./dx)./sqrt(1+(draft./dx)^2)))...
            +trapz(x(indl:indr)-Xc,(x(indl:indr)-X0-Xc).*((dxPsiP_x(indl:indr).^2+dzPsiP_x.^2)/2-dzPsiP_x.*nu.z(1,indl:indr))));
    end
    
    if strcmpi(Idmotion,'Surge')|| strcmpi(Idmotion,'Free')
        psiS=psiS.';psiS(isnan(psiS))=0;
        FpsiEvS(indl)=FevLS;FpsiEvS(indr)=FevRS;
        
        FpsiS=tau1S.*(D-draft);
        FpsiS=FpsiS.*chiSWL;
    end
    if strcmpi(Idmotion,'Heave')|| strcmpi(Idmotion,'Free')
        psiH=psiH.';psiH(isnan(psiH))=0;
        FpsiEvH(indl)=FevLH;FpsiEvH(indr)=FevRH;
        FpsiH=-(x-Xc)+tau1H*(D-draft)+(dxD./(D-draft)).*((D.^2+draft^2)/2-draft.*D+(-(D-draft).^2/3+(x-Xc).^2)./2);
        FpsiH=FpsiH.*chiSWL;
    end
    if strcmpi(Idmotion,'Pitch')|| strcmpi(Idmotion,'Free')
        psiP=psiP.';psiP(isnan(psiP))=0;
        FpsiEvP(indl)=FevLP;FpsiEvP(indr)=FevRP;
        
        FpsiP=(D-draft).^2/6-(x-Xc).^2/2+tau1P*(D-draft)+(dxD./(D-draft)).*((x-Xc).*((D.^2+draft^2)/2-draft.*D)+(-(D-draft).^2/3.*(x-Xc)+(x-Xc).^3/3)./2);
        FpsiP=FpsiP.*chiSWL;
    end
    Ma.z=Ma33;Ma.x=Ma11;
    Ma.theta=Ma55;
    if strcmpi(Idmotion,'Free')
        Ma.S(1).Mat=[Ma11 Ma13 Ma15;Ma31 Ma33 Ma35;Ma51 Ma53 Ma55];
        dampC.S(1).Mat=[B11 B13 B15;B31 B33 B35;B51 B53 B55];
        %Ma.Mat=diag([Ma11 Ma33 Ma55]);
    else
        Ma.S(1).Mat=diag([Ma11 Ma33 Ma55]);
        dampC.S(1).Mat=diag([B11 B33 B55]);
    end
    Fpsi.surge=FpsiS;Fpsi.heave=FpsiH;
    Fpsi.pitch=FpsiP;
    FpsiEv.surge=FpsiEvS;FpsiEv.heave=FpsiEvH;
    FpsiEv.pitch=FpsiEvP;
    psi.surge=psiS;psi.heave=psiH;
    psi.pitch=psiP;
    Lpsiev.surge=LpsievS;Lpsiev.heave=LpsievH;
    Lpsiev.pitch=LpsievP;
    
    
