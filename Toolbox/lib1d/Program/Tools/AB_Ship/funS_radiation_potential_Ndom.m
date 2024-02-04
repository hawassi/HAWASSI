function [psi,Ma,dampC,Lpsiev,Fpsi,FpsiEv,dxi_Krad]=funS_radiation_potential_Ndom(shippar,x,xShip,Sdraft,sZ,chiXcZ0,chiWL,nu,Smotion,D,Nship,etaL,etaR,dxetaL,dxetaR)

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
    Evmodes=shippar.Evmode;
    b=shippar.form.Slength(1)/2;
    kappaS=shippar.form.kappa;
    kappaS(1,:)=-kappaS(1,:);
    
    draft=Sdraft(ii)-sZ(ii);
   
    Zc=shippar.form.XZrot0(ii,2);%% for pitch
    
    chiship    =shippar.form.chi;
    chiwaterline=shippar.form.chiWl;
    chiWl      =chiwaterline(:,end).';
    chiS      =chiship(:,end).';
    chiSWl     =chiS+chiWl;

    etaWl=[etaL etaR];
    dxetaWl=[dxetaL dxetaR];
    
    
    shipformOnGrid=shippar.form.shapeXcZc;
    taug=shippar.tangenflowGlobal.';
    dxtaug=gradient(taug,dx);
    
        Ndom=shippar.form.shapedat(1).Ndom;
        
        XshipL=shippar.form.shapedat(1).XshipL;
        ZshipL=shippar.form.shapedat(1).ZshipL;
        XshipR=shippar.form.shapedat(1).XshipR;
        ZshipR=shippar.form.shapedat(1).ZshipR;
        dxZshipL=shippar.form.shapedat(1).dxZshipL;
        dxZshipR=shippar.form.shapedat(1).dxZshipR;
        
        Xbdy=shippar.form.shapedat(1).Xbdy;
        dr=shippar.form.shapedat(1).dr_approx-sZ(ii);

        indxBdy=shippar.form.shapedat(1).indxBdy;
         indx1=indxBdy(1);indx2=indxBdy(end);
         shipform=shipformOnGrid(end,:);
        
          
        Dbdy=D(indxBdy);
        dxDbdy=dxD(indxBdy);
        tauGbdy=taug(indxBdy);
        dxtauGbdy=dxtaug(indxBdy);
    
         DZeta=D+shipform.';
        %dxDZeta=funOprt_FDGradient1d(DZeta,dx,2);
        DZetabdy=DZeta(indxBdy);%[DZeta(indx1) DZeta(indxlc) DZeta(indc) DZeta(indxcr) DZeta(indx2)];
        dxDZetabdy=zeros(size(DZetabdy));%%dxDZeta(indxBdy);%[dxDbdy(1) dxDZeta(indxlc) dxDZeta(indc) dxDZeta(indxcr) dxDbdy(5)];

     [Amat,BvectS,BvectH,BvectP] =fun_construct_matrices_radiation_Ndom(Ndom,...
                0,0,tauGbdy,dxtauGbdy,Dbdy,dxDbdy,DZetabdy,dxDZetabdy,dr,...
                Evmodes,kappaS,Idmotion,Xc,Zc,etaWl,dxetaWl,Xbdy);

                
  
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
    GammaH=zeros(size(x));
    GammaP=zeros(size(x));
    
    
     %%%for calculating d\xi K(\phi,\xi)
     DZL=DL+ZshipL;
     DZR=DR+ZshipR;
%      dzQLz0=-kappa0L.*sin(kappa0L.*DZL)./cos(kappa0L.*HL);
%      dzQRz0=-kappa0R.*sin(kappa0R.*DZR)./cos(kappa0R.*HR);
%   
    InQ0TEtaL=fun_IntegQn(kappaS(1,1),etaL,Dbdy(1),-dr(1),etaL);
    InQ0TEtaR=fun_IntegQn(kappaS(1,2),etaR,Dbdy(end),-dr(end),etaR);

    
    if strcmpi(Idmotion,'Surge')|| strcmpi(Idmotion,'Free')
        CoefS=Amat\BvectS;
        psiS=zeros(size(chiS));
        Nm=Ndom*2-2;
        idAd=0;idCoef1=3;idC=Ndom/2;idAdE=0;idCoef2=Nm-1;
        f11d=CoefS(1)*InQ0TEtaL-CoefS(Nm)*InQ0TEtaR;
        
        InZQ0TEtaL=fun_IntegZQn(kappaS(1,1),etaL,Dbdy(1),-dr(1),etaL);
        InZQ0TEtaR=fun_IntegZQn(kappaS(1,2),etaR,Dbdy(end),-dr(end),etaR);
        f151d=-CoefS(1)*(InZQ0TEtaL-Zc*InQ0TEtaL)+CoefS(Nm)*(InZQ0TEtaR-Zc*InQ0TEtaR);
        
        for jj=1:(Ndom-2)/2
            if jj==(Ndom-2)/2, idAdE=1; end
            int1=indxBdy(jj)+idAd:indxBdy(jj+1);
            psiS(int1)=psiS(int1)+CoefS(idCoef1).*taug(int1)+CoefS(idCoef1-1);
            int2=indxBdy(end-jj)+idAdE:indxBdy(end-jj+1)-idAd;
            psiS(int2)=psiS(int2)+CoefS(idCoef2).*taug(int2)+CoefS(idCoef2-1);
            
            if jj<(Ndom-2)/2
                Intvarphi0Zetad_dprev=fun_IntegVarphin(0,Dbdy(jj+1),-dr(jj+1),-dr(jj));
                f11d=f11d+(CoefS(idCoef1).*taug(indxBdy(jj+1))+CoefS(idCoef1-1))*Intvarphi0Zetad_dprev;
                
                Intvarphi0Zetad_dnext=fun_IntegVarphin(0,Dbdy(end-jj),-dr(end-jj),-dr(end-jj+1));
                f11d=f11d-(CoefS(idCoef2).*taug(indxBdy(end-jj))+CoefS(idCoef2-1))*Intvarphi0Zetad_dnext;
                
                IntZvarphi0Zetad_dprev=fun_IntegZVarphin(0,Dbdy(jj+1),-dr(jj+1),-dr(jj));
                IntZvarphi0Zetad_dnext=fun_IntegZVarphin(0,Dbdy(end-jj),-dr(end-jj),-dr(end-jj+1));
                
                f151d=f151d-(CoefS(idCoef1).*taug(indxBdy(jj+1))+CoefS(idCoef1-1))*(IntZvarphi0Zetad_dprev-Zc*Intvarphi0Zetad_dprev);
                f151d=f151d+(CoefS(idCoef2).*taug(indxBdy(end-jj))+CoefS(idCoef2-1))*(IntZvarphi0Zetad_dnext-Zc*Intvarphi0Zetad_dnext);
            end
            idAd=1;idCoef1=idCoef1+2; idCoef2=idCoef2-2;
        end
        
        psiS=psiS.*chiS;
        
        psiS(indxBdy(1))=CoefS(1);
        psiS(indxBdy(end))=CoefS(Nm);
        
%         %%%
%         dPhi1Sdz_zL=real(ALS.*dzQLz0);dPhi1Sdx_zL=1;
%         dPhi3Sdz_zR=real(ARS.*dzQRz0);dPhi3Sdx_zR=1;
%         %%%
    end
    
    if strcmpi(Idmotion,'Heave')|| strcmpi(Idmotion,'Free')
        CoefH=Amat\BvectH;
        GammaH(indx1:indx2)=((D(indx1:indx2)+shipform(indx1:indx2)).^2-(x(indx1:indx2)-Xc).^2)./2./(D(indx1:indx2)+shipform(indx1:indx2));
           
      
        psiH=GammaH.*chiSWl;
        psiH(indxBdy(1))=GammaH(indxBdy(1)+1);
        psiH(indxBdy(end))=GammaH(indxBdy(end)-1);
        Nm=Ndom*2-2;
        idAd=0;idCoef1=3;idC=Ndom/2;idAdE=0;idCoef2=Nm-1;
    
        
        f31d=CoefH(1)*InQ0TEtaL-CoefH(Nm)*InQ0TEtaR;
        
        InZQ0TEtaL=fun_IntegZQn(kappaS(1,1),etaL,Dbdy(1),-dr(1),etaL);
        InZQ0TEtaR=fun_IntegZQn(kappaS(1,2),etaR,Dbdy(end),-dr(end),etaR);
        
        f351d=-CoefH(1)*(InZQ0TEtaL-Zc*InQ0TEtaL)+CoefH(Nm)*(InZQ0TEtaR-Zc*InQ0TEtaR);
        
        
        for jj=1:(Ndom-2)/2
            if jj==(Ndom-2)/2, idAdE=1; end
            int1=indxBdy(jj)+idAd:indxBdy(jj+1);
            psiH(int1)=psiH(int1)+CoefH(idCoef1).*taug(int1)+CoefH(idCoef1-1);
            int2=indxBdy(end-jj)+idAdE:indxBdy(end-jj+1)-idAd;
            psiH(int2)=psiH(int2)+CoefH(idCoef2).*taug(int2)+CoefH(idCoef2-1);
            
            if jj<(Ndom-2)/2
                Intvarphi0Zetad_dprev=fun_IntegVarphin(0,Dbdy(jj+1),-dr(jj+1),-dr(jj));
                f31d=f31d+(CoefH(idCoef1).*taug(indxBdy(jj+1))+CoefH(idCoef1-1))*Intvarphi0Zetad_dprev;
                
                Intvarphi0Zetad_dnext=fun_IntegVarphin(0,Dbdy(end-jj),-dr(end-jj),-dr(end-jj+1));
                f31d=f31d-(CoefH(idCoef2).*taug(indxBdy(end-jj))+CoefH(idCoef2-1))*Intvarphi0Zetad_dnext;
                
                IntZvarphi0Zetad_dprev=fun_IntegZVarphin(0,Dbdy(jj+1),-dr(jj+1),-dr(jj));
                IntZvarphi0Zetad_dnext=fun_IntegZVarphin(0,Dbdy(end-jj),-dr(end-jj),-dr(end-jj+1));
                
                f351d=f351d-(CoefH(idCoef1).*taug(indxBdy(jj+1))+CoefH(idCoef1-1))*(IntZvarphi0Zetad_dprev-Zc*Intvarphi0Zetad_dprev);
                f351d=f351d+(CoefH(idCoef2).*taug(indxBdy(end-jj))+CoefH(idCoef2-1))*(IntZvarphi0Zetad_dnext-Zc*Intvarphi0Zetad_dnext);
            end
            
            
            idAd=1;idCoef1=idCoef1+2; idCoef2=idCoef2-2;
        end
        psiH=psiH.*chiSWl;
        
%          %%%for calculating d\xi K(\phi,\xi)
%         dxGamma3_x=-(x-Xc)./(D-draft)+dxD.*(1-((D-draft)-(x-Xc).^2)./2./(D-draft));
%         dxPsiH_x=real(dxGamma3_x+tau1H);
%         dzPsiH_x=1 ;
%         %%%%%%%%%%%
    end
    
    
    if strcmpi(Idmotion,'Pitch')|| strcmpi(Idmotion,'Free')
        CoefP=Amat\BvectP;
        
        GammaP(indx1:indx2)=((D(indx1:indx2)+shipform(indx1:indx2)).^2.*(x(indx1:indx2)-Xc)-(x(indx1:indx2)-Xc).^3/3)./2./(D(indx1:indx2)+shipform(indx1:indx2));
               
        psiP=GammaP.*chiSWL;
        psiP(indx1)=GammaP(indx1+1);
        psiP(indx2)=GammaP(indx2-1);
        %
        Nm=Ndom*2-2;
        idAd=0;idCoef1=3;idC=Ndom/2;idAdE=0;idCoef2=Nm-1;
        InZQ0TEtaL=fun_IntegZQn(kappaS(1,1),etaL,Dbdy(1),-dr(1),etaL);
        InZQ0TEtaR=fun_IntegZQn(kappaS(1,2),etaR,Dbdy(end),-dr(end),etaR);
        f51d=CoefP(1)*InQ0TEtaL-CoefP(Nm)*InQ0TEtaR;
        
        f551d=-CoefP(1)*(InZQ0TEtaL-Zc*InQ0TEtaL)+CoefP(Nm)*(InZQ0TEtaR-Zc*InQ0TEtaR);
        
        for jj=1:(Ndom-2)/2
            if jj==(Ndom-2)/2, idAdE=1; end
            int1=indxBdy(jj)+idAd:indxBdy(jj+1);
            psiP(int1)=psiP(int1)+CoefP(idCoef1).*taug(int1)+CoefP(idCoef1-1);
            int2=indxBdy(end-jj)+idAdE:indxBdy(end-jj+1)-idAd;
            psiP(int2)=psiP(int2)+CoefP(idCoef2).*taug(int2)+CoefP(idCoef2-1);
            
            
            if jj<(Ndom-2)/2
                Intvarphi0Zetad_dprev=fun_IntegVarphin(0,Dbdy(jj+1),-dr(jj+1),-dr(jj));
                IntZvarphi0Zetad_dprev=fun_IntegZVarphin(0,Dbdy(jj+1),-dr(jj+1),-dr(jj));
                IntGammaPZZ0_dprev=fun_IntGammaPZZ0(Dbdy(jj+1),Xbdy(jj+1)-Xc,Zc,-dr(jj+1),-dr(jj));
               
                f51d=f51d+(CoefP(idCoef1).*taug(indxBdy(jj+1))+CoefP(idCoef1-1))*Intvarphi0Zetad_dprev;
                f551d=f551d-IntGammaPZZ0_dprev-(CoefP(idCoef1).*taug(indxBdy(jj+1))+CoefP(idCoef1-1))*(IntZvarphi0Zetad_dprev-Zc*Intvarphi0Zetad_dprev);
                
                Intvarphi0Zetad_dnext=fun_IntegVarphin(0,Dbdy(end-jj),-dr(end-jj),-dr(end-jj+1));
                IntZvarphi0Zetad_dnext=fun_IntegZVarphin(0,Dbdy(end-jj),-dr(end-jj),-dr(end-jj+1));
                IntGammaPZZ0_dnext=fun_IntGammaPZZ0(Dbdy(end-jj),Xbdy(end-jj)-Xc,Zc,-dr(end-jj),-dr(end-jj+1));
                f51d=f51d-(CoefP(idCoef2).*taug(indxBdy(end-jj))+CoefP(idCoef2-1))*Intvarphi0Zetad_dnext;
                f551d=f551d+IntGammaPZZ0_dnext+(CoefP(idCoef2).*taug(indxBdy(end-jj))+CoefP(idCoef2-1))*(IntZvarphi0Zetad_dnext-Zc*Intvarphi0Zetad_dnext);
                
               
            end
            
            idAd=1;idCoef1=idCoef1+2; idCoef2=idCoef2-2;
        end
        psiP=psiP.*chiSWL;
        
%          %%%for calculating d\xi K(\phi,\xi)
%         dxGamma5_x=((D-draft).^2-(x-Xc).^2)./2./(D-draft)+dxD.*((x-Xc)-((D-draft).*(x-Xc)-(x-Xc).^3/3)./2./(D-draft));
%         dxPsiP_x=real(dxGamma5_x+tau1P);
%         dzPsiP_x=1 ;
%         dPhi1Pdz_zL=real(ALP.*dzQLz0);dPhi1Pdx_zL=1;
%         dPhi3Pdz_zR=real(ARP.*dzQRz0);dPhi3Pdx_zR=1;
%         
    end
    
    
    
    
%     if strcmpi(Idmotion,'Free')    
%           f13= 2*tau0S*b;
%           f151=ARS.*(IZQ0TetaR-Zc.*IntQ0TetaR)-ALS.*(IZQ0TetaL-Zc.*IntQ0TetaL);%
%           f152=(-tau0S*(Xc*2*b)+tau1S*2*b^3/3);
%           f31=ALH.*IntQ0TetaL-ARH*IntQ0TetaR;   
%           f351=ARH.*(IZQ0TetaR-Zc.*IntQ0TetaR)-ALH.*(IZQ0TetaL-Zc.*IntQ0TetaL);%;
%           f352=trapz(x(indl:indr)-Xc,(x(indl:indr)-Xc).*psiH(indl:indr));
%           f51=ALP.*IntQ0TetaL-ARP*IntQ0TetaR;
%           f53= trapz(x(indl:indr)-Xc,psiP(indl:indr));
%     end
    
    if Evmodes>0
        for m=1:Evmodes
            mu_m=m*pi./(DZeta);
            
             InQmTEtaL=fun_IntegQn(kappaS(m+1,1),etaL,Dbdy(1),-dr(1),etaL);
             InQmTEtaR=fun_IntegQn(kappaS(m+1,2),etaR,Dbdy(end),-dr(end),etaR);
           
            if strcmpi(Idmotion,'Surge')|| strcmpi(Idmotion,'Free')
                psiSev=zeros(size(x));
                %
                idAd=0;idAdE=0; idCoef1=Nm*m+2; idCoef2=Nm*m+Nm-2;
                f11d=f11d+CoefS(Nm*m+1)*InQmTEtaL-CoefS(Nm*(m+1))*InQmTEtaR;
                
                InZQmTEtaL=fun_IntegZQn(kappaS(m+1,1),etaL,Dbdy(1),-dr(1),etaL);
                InZQmTEtaR=fun_IntegZQn(kappaS(m+1,2),etaR,Dbdy(end),-dr(end),etaR);
                f151d=f151d-CoefS(Nm*m+1)*(InZQmTEtaL-Zc*InQmTEtaL)+CoefS(Nm*(m+1))*(InZQmTEtaR-Zc*InQmTEtaR);
                
                for jj=1:(Ndom-2)/2
                    if jj==(Ndom-2)/2, idAdE=1; end
                    int1=indxBdy(jj)+idAd:indxBdy(jj+1);
                    psiSev(int1)=((CoefS(idCoef1).*exp(mu_m(int1).*(x(int1)-Xbdy(jj)))...
                        +CoefS(idCoef1+1).*exp(-mu_m(int1).*(x(int1)-Xbdy(jj+1)))).*cos(mu_m(int1).*(DZeta(int1))));
                    int2=indxBdy(end-jj)+idAdE:indxBdy(end-jj+1)-idAd;
                    psiSev(int2)=((CoefS(idCoef2).*exp(mu_m(int2).*(x(int2)-Xbdy(end-jj)))...
                        +CoefS(idCoef2+1).*exp(-mu_m(int2).*(x(int2)-Xbdy(end-jj+1)))).*cos(mu_m(int2).*(DZeta(int2))));
                    
                    if jj<(Ndom-2)/2
                        mu_mLd=m*pi./(Dbdy(jj)-dr(jj));
                        mu_mRd=m*pi./(Dbdy(end-jj+1)-dr(end-jj+1));
                        IntvarphimZetad_dprev=fun_IntegVarphin(mu_mLd,Dbdy(jj+1),-dr(jj+1),-dr(jj));
                        f11d=f11d+(CoefS(idCoef1).*exp(mu_mLd.*(Xbdy(jj+1)-Xbdy(jj)))...
                            +CoefS(idCoef1+1)).*IntvarphimZetad_dprev;
                        IntvarphimZetad_dnext=fun_IntegVarphin(mu_mRd,Dbdy(end-jj),-dr(end-jj),-dr(end-jj+1));
                        f11d=f11d-(CoefS(idCoef2)+CoefS(idCoef2+1).*exp(-mu_mRd.*(Xbdy(end-jj)-Xbdy(end-jj+1))))*IntvarphimZetad_dnext;
                        
                        IntZvarphimZetad_dprev=fun_IntegZVarphin(mu_mLd,Dbdy(jj+1),-dr(jj+1),-dr(jj));
                        f151d=f151d-(CoefS(idCoef1).*exp(mu_mLd.*(Xbdy(jj+1)-Xbdy(jj)))...
                            +CoefS(idCoef1+1)).*(IntZvarphimZetad_dprev-Zc*IntvarphimZetad_dprev);
                        %
                        IntZvarphimZetad_dnext=fun_IntegZVarphin(mu_mRd,Dbdy(end-jj),-dr(end-jj),-dr(end-jj+1));
                        
                        f151d=f151d+(CoefS(idCoef2)+CoefS(idCoef2+1).*...
                            exp(-mu_mRd.*(Xbdy(end-jj)-Xbdy(end-jj+1))))*(IntZvarphimZetad_dnext-Zc*IntvarphimZetad_dnext);
                        %
                    end
                    
                    idAd=1;idCoef1=idCoef1+2;idCoef2=idCoef2-2;
                end
                psiSev=psiSev.*chiS;
                psiSev(indxBdy(1))=CoefS(Nm*m+1);
                psiSev(indxBdy(end))=CoefS(Nm*(m+1));
                psiS=psiS+psiSev;
                
                %%%%
%                 LpsievS(1:indl)=LpsievS(1:indl)-epsS1m*exp(kappaL(m).*(x(1:indl)-Xc+b)).*kappaL(m).*tan(kappaL(m).*D(1:indl));
%                 LpsievS(indr:end)=LpsievS(indr:end)-epsS3m*exp(-kappaR(m).*(x(indr:end)-Xc-b)).*kappaR(m).*tan(kappaR(m).*D(indr:end));
%                 FevLS=FevLS+epsS1m*tan(kappaL(m).*HL);
%                 FevRS=FevRS-epsS3m*tan(kappaR(m).*HR);
                
                %%%for calculating d\xi K(\phi,\xi)
%                 dzQLzm=-kappaL(m).*sin(kappaL(m).*DZL)./cos(kappaL(m).*HL);
%                 dzQRzm=-kappaR(m).*sin(kappaR(m).*DZR)./cos(kappaR(m).*HR);
%    
%                 dPhi1Sdz_zL=dPhi1Sdz_zL+epsS1m.*dzQLzm;
%                 dPhi3Sdz_zR=dPhi3Sdz_zR+epsS3m.*dzQRzm;
            end
            
            if strcmpi(Idmotion,'Heave')|| strcmpi(Idmotion,'Free')
                psiHev=zeros(size(x));
                %
                idAd=0;idAdE=0; idCoef1=Nm*m+2; idCoef2=Nm*m+Nm-2;
               
                f31d=f31d+CoefH(Nm*m+1)*InQmTEtaL-CoefH(Nm*(m+1))*InQmTEtaR;
                
                InZQmTEtaL=fun_IntegZQn(kappaS(m+1,1),etaL,Dbdy(1),-dr(1),etaL);
                InZQmTEtaR=fun_IntegZQn(kappaS(m+1,2),etaR,Dbdy(end),-dr(end),etaR);
                f351d=f351d-CoefH(Nm*m+1)*(InZQmTEtaL-Zc*InQmTEtaL)+CoefH(Nm*(m+1))*(InZQmTEtaR-Zc*InQmTEtaR);

                
                for jj=1:(Ndom-2)/2
                    if jj==(Ndom-2)/2, idAdE=1; end
                    int1=indxBdy(jj)+idAd:indxBdy(jj+1);
                    psiHev(int1)=((CoefH(idCoef1).*exp(mu_m(int1).*(x(int1)-Xbdy(jj)))...
                        +CoefH(idCoef1+1).*exp(-mu_m(int1).*(x(int1)-Xbdy(jj+1)))).*cos(mu_m(int1).*(DZeta(int1))));
                    int2=indxBdy(end-jj)+idAdE:indxBdy(end-jj+1)-idAd;
                    psiHev(int2)=((CoefH(idCoef2).*exp(mu_m(int2).*(x(int2)-Xbdy(end-jj)))...
                        +CoefH(idCoef2+1).*exp(-mu_m(int2).*(x(int2)-Xbdy(end-jj+1)))).*cos(mu_m(int2).*(DZeta(int2))));
                    if jj<(Ndom-2)/2
                        mu_mLd=m*pi./(Dbdy(jj)-dr(jj));
                        mu_mRd=m*pi./(Dbdy(end-jj+1)-dr(end-jj+1));
                        IntvarphimZetad_dprev=fun_IntegVarphin(mu_mLd,Dbdy(jj+1),-dr(jj+1),-dr(jj));
                        f31d=f31d+(CoefH(idCoef1).*exp(mu_mLd.*(Xbdy(jj+1)-Xbdy(jj)))...
                            +CoefH(idCoef1+1)).*IntvarphimZetad_dprev;
                        IntvarphimZetad_dnext=fun_IntegVarphin(mu_mRd,Dbdy(end-jj),-dr(end-jj),-dr(end-jj+1));
                        f31d=f31d-(CoefH(idCoef2)+CoefH(idCoef2+1).*exp(-mu_mRd.*(Xbdy(end-jj)-Xbdy(end-jj+1))))*IntvarphimZetad_dnext;
                        
                        IntZvarphimZetad_dprev=fun_IntegZVarphin(mu_mLd,Dbdy(jj+1),-dr(jj+1),-dr(jj));
                        f351d=f351d-(CoefH(idCoef1).*exp(mu_mLd.*(Xbdy(jj+1)-Xbdy(jj)))...
                            +CoefH(idCoef1+1)).*(IntZvarphimZetad_dprev-Zc*IntvarphimZetad_dprev);
                        %
                        IntZvarphimZetad_dnext=fun_IntegZVarphin(mu_mRd,Dbdy(end-jj),-dr(end-jj),-dr(end-jj+1));
                        
                        f351d=f351d+(CoefH(idCoef2)+CoefH(idCoef2+1).*...
                            exp(-mu_mRd.*(Xbdy(end-jj)-Xbdy(end-jj+1))))*(IntZvarphimZetad_dnext-Zc*IntvarphimZetad_dnext);
                        %
                    end
     
                    idAd=1;idCoef1=idCoef1+2;idCoef2=idCoef2-2;
                end
                psiHev=psiHev.*chiSWL;
                psiH=psiH+psiHev;
                
                
                %%%%%%%%%%%%
%                 LpsievH(1:indl)=LpsievH(1:indl)-epsH1m*exp(kappaL(m).*(x(1:indl)-Xc+b)).*kappaL(m).*tan(kappaL(m).*D(1:indl));
%                 LpsievH(indr:end)=LpsievH(indr:end)-epsH3m*exp(-kappaR(m).*(x(indr:end)-Xc-b))*kappaR(m).*tan(kappaR(m).*D(indr:end));
%                 FevLH=FevLH+epsH1m*tan(kappaL(m)*HL);
%                 FevRH=FevRH-epsH3m*tan(kappaR(m)*HR);
                
%                 dxPsiHev=mu_m.*(1-(x-Xc).*dxD./(D-draft)).*(epsH21m.*exp(mu_m.*((x-Xc)+b))-epsH22m.*exp(-mu_m.*((x-Xc)-b))).*cos(mu_m.*(D-draft));
%                 dxPsiH_x=dxPsiH_x+real(dxPsiHev);
            end
           
            if strcmpi(Idmotion,'Pitch')|| strcmpi(Idmotion,'Free')
                psiPev=zeros(size(x));
                
                InQmTEtaL=fun_IntegQn(kappaS(m+1,1),etaL,Dbdy(1),-dr(1),etaL);
                InQmTEtaR=fun_IntegQn(kappaS(m+1,2),etaR,Dbdy(end),-dr(end),etaR);
                InZQmTEtaL=fun_IntegZQn(kappaS(m+1,1),etaL,Dbdy(1),-dr(1),etaL);
                InZQmTEtaR=fun_IntegZQn(kappaS(m+1,2),etaR,Dbdy(end),-dr(end),etaR);
                
                f51d=f51d+CoefP(Nm*m+1)*InQmTEtaL-CoefP(Nm*(m+1))*InQmTEtaR;
                f551d=f551d-CoefP(Nm*m+1)*(InZQmTEtaL-Zc*InQmTEtaL)+CoefP(Nm*(m+1))*(InZQmTEtaR-Zc*InQmTEtaR);
                
                idAd=0;idAdE=0; idCoef1=Nm*m+2; idCoef2=Nm*m+Nm-2;
                for jj=1:(Ndom-2)/2
                    if jj==(Ndom-2)/2, idAdE=1; end
                    int1=indxBdy(jj)+idAd:indxBdy(jj+1);
                    psiPev(int1)=((CoefP(idCoef1).*exp(mu_m(int1).*(x(int1)-Xbdy(jj)))...
                        +CoefP(idCoef1+1).*exp(-mu_m(int1).*(x(int1)-Xbdy(jj+1)))).*cos(mu_m(int1).*(DZeta(int1))));
                    int2=indxBdy(end-jj)+idAdE:indxBdy(end-jj+1)-idAd;
                    psiPev(int2)=((CoefP(idCoef2).*exp(mu_m(int2).*(x(int2)-Xbdy(end-jj)))...
                        +CoefP(idCoef2+1).*exp(-mu_m(int2).*(x(int2)-Xbdy(end-jj+1)))).*cos(mu_m(int2).*(DZeta(int2))));
                    
                    
                    if jj<(Ndom-2)/2
                        mu_mLd=m*pi./(Dbdy(jj)-dr(jj));
                        mu_mRd=m*pi./(Dbdy(end-jj+1)-dr(end-jj+1));
                        IntvarphimZetad_dprev=fun_IntegVarphin(mu_mLd,Dbdy(jj+1),-dr(jj+1),-dr(jj));
                        IntZvarphimZetad_dprev=fun_IntegZVarphin(mu_mLd,Dbdy(jj+1),-dr(jj+1),-dr(jj));
                        f551d=f551d-(CoefP(idCoef1).*exp(mu_mLd.*(Xbdy(jj+1)-Xbdy(jj)))...
                            +CoefP(idCoef1+1)).*(IntZvarphimZetad_dprev-Zc*IntvarphimZetad_dprev);
                        
                        f51d=f51d+(CoefP(idCoef1).*exp(mu_mLd.*(Xbdy(jj+1)-Xbdy(jj)))...
                            +CoefP(idCoef1+1)).*IntvarphimZetad_dprev;
                        IntvarphimZetad_dnext=fun_IntegVarphin(mu_mRd,Dbdy(end-jj),-dr(end-jj),-dr(end-jj+1));
                        IntZvarphimZetad_dnext=fun_IntegZVarphin(mu_mRd,Dbdy(end-jj),-dr(end-jj),-dr(end-jj+1));
                        
                        f551d=f551d+(CoefP(idCoef2)+CoefP(idCoef2+1).*...
                            exp(-mu_mRd.*(Xbdy(end-jj)-Xbdy(end-jj+1))))*(IntZvarphimZetad_dnext-Zc*IntvarphimZetad_dnext);
                   
                        f51d=f51d-(CoefP(idCoef2)+CoefP(idCoef2+1).*exp(-mu_mRd.*(Xbdy(end-jj)-Xbdy(end-jj+1))))*IntvarphimZetad_dnext;
                    end
                    idAd=1;idCoef1=idCoef1+2;idCoef2=idCoef2-2;
                end
                psiPev=psiPev.*chiSWL;
                psiP=psiP+psiPev;
                
                
%                 LpsievP(1:indl)=LpsievP(1:indl)-epsP1m*exp(kappaL(m).*(x(1:indl)-Xc+b)).*kappaL(m).*tan(kappaL(m).*D(1:indl));
%                 LpsievP(indr:end)=LpsievP(indr:end)-epsP3m.*exp(-kappaR(m).*(x(indr:end)-Xc-b)).*kappaR(m).*tan(kappaR(m).*D(indr:end));
%                 FevLP=FevLP+epsP1m*tan(kappaL(m)*HL);
%                 FevRP=FevRP-epsP3m*tan(kappaR(m)*HR);
                
                
                %%%for calculating d\xi K(\phi,\xi)
%                 dzQLzm=-kappaL(m).*sin(kappaL(m).*DZL)./cos(kappaL(m).*HL);
%                 dzQRzm=-kappaR(m).*sin(kappaR(m).*DZR)./cos(kappaR(m).*HR);
%                 
%                 dPhi1Pdz_zL=dPhi1Pdz_zL+epsP1m.*dzQLzm;
%                 dPhi3Pdz_zR=dPhi3Pdz_zR+epsP3m.*dzQRzm;
%                 
%                 dxPsiPev=mu_m.*(1-(x-Xc).*dxD./(D-draft)).*(epsP21m.*exp(mu_m.*((x-Xc)+b))-epsP22m.*exp(-mu_m.*((x-Xc)-b))).*cos(mu_m.*(D-draft));
%                 dxPsiP_x=dxPsiP_x+real(dxPsiPev);
            end
            
            
            
            
%             if strcmpi(Idmotion,'Free')
%                  f13=f13+trapz(x(indl:indr)-Xc,psiSev(indl:indr));
%                  f151=f151+epsS3m.*(IZQmTetaR-Zc.*IntQmTetaR)-epsS1m.*(IZQmTetaL-Zc.*IntQmTetaL);
%                  f152=f152+trapz(x(indl:indr)-Xc,(x(indl:indr)-Xc).*psiSev(indl:indr));
%                  f31=f31+(epsH1m.*IntQmTetaL-epsH3m.*IntQmTetaR);
%                  f351=f351+epsH3m.*(IZQmTetaR-Zc.*IntQmTetaR)-epsH1m.*(IZQmTetaL-Zc.*IntQmTetaL);%
%                  f352=f352+trapz(x(indl:indr)-Xc,(x(indl:indr)-Xc).*psiHev(indl:indr));
%                  f51=f51+(epsP1m.*IntQmTetaL-epsP3m.*IntQmTetaR);
%                  f53=f53+trapz(x(indl:indr)-Xc,psiPev(indl:indr));
%             end
        end
    end
    %      figure;
    %      plot(x,psi);
    dxi_Krad.z(1)      = 0;
    dxi_Krad.x(1)      = 0;
    dxi_Krad.theta(1)   = 0;
    
    if strcmpi(Idmotion,'Surge')
        Ma33=0;Ma55=0;
        Ma11=real(f11d);
        B11=imag(f11d);B33=0;B55=0;
%         dxi_Krad.x(ii)= real(trapz(ZshipL,(dPhi1Sdx_zL.^2+dPhi1Sdz_zL.^2)./2-dPhi1Sdz_zL.*(dPhi1Sdx_zL.*draft./dx)./sqrt(1+(draft./dx)^2)) ...
%             -trapz(ZshipR,(dPhi3Sdx_zR.^2+dPhi3Sdz_zR.^2)./2+dPhi3Sdz_zR.*(dPhi3Sdx_zR.*draft./dx)./sqrt(1+(draft./dx)^2)));
    elseif strcmpi(Idmotion,'Heave')
        f33=trapz(x(indxBdy(1):indxBdy(end)),psiH(indxBdy(1):indxBdy(end)));
              
        Ma33=real(f33);B33=imag(f33);
        Ma11=0;Ma55=0;B55=0;B11=0;
        % zMa=trapz(x(indl:indr),psi(indl:indr).*nutild.z(1,indl:indr));
        %      [zMa zMaN]
        
%         dxi_Krad.z(1)     =  real(trapz(x(indl:indr)-Xc,(dxPsiH_x(indl:indr).^2+dzPsiH_x.^2)/2-dzPsiH_x.*nu.z(1,indl:indr)));%trapz(x(indl:indr)-Xc,(dxPsiH_x(indl:indr).^2-dzPsiH_x.^2)/2);%
    elseif strcmpi(Idmotion,'Pitch')
        f552=trapz(x(indxBdy(1):indxBdy(end)),(x(indxBdy(1):indxBdy(end))-Xc).*psiP(indxBdy(1):indxBdy(end)));
        f55=f551d+f552;
        Ma55=real(f55);B55=imag(f55);B33=0;B11=0;
        %     thetaMaN=trapz(x(indl:indr),psi(indl:indr).*nutild.theta(1,indl:indr))
        Ma11=0;Ma33=0;
%         dxi_Krad.theta(1)  = real(-trapz(ZshipL,(ZshipL-Zc).*((dPhi1Pdx_zL.^2+dPhi1Pdz_zL.^2)./2-dPhi1Pdz_zL.*(dPhi1Pdx_zL.*draft./dx)./sqrt(1+(draft./dx)^2))) ...
%             +trapz(ZshipR,(ZshipL-Zc).*((dPhi3Pdx_zR.^2+dPhi3Pdz_zR.^2)./2+dPhi3Pdz_zR.*(dPhi3Pdx_zR.*draft./dx)./sqrt(1+(draft./dx)^2)))...
%             +trapz(x(indl:indr)-Xc,(x(indl:indr)-Xc-Xc).*((dxPsiP_x(indl:indr).^2+dzPsiP_x.^2)/2-dzPsiP_x.*nu.z(1,indl:indr))));      
    elseif strcmpi(Idmotion,'Free')
        f13= trapz(x(indxBdy(1):indxBdy(end)),psiS(indxBdy(1):indxBdy(end)));
        f152=trapz(x(indxBdy(1):indxBdy(end)),(x(indxBdy(1):indxBdy(end))-Xc).*psiS(indxBdy(1):indxBdy(end)));
        f552=trapz(x(indxBdy(1):indxBdy(end)),(x(indxBdy(1):indxBdy(end))-Xc).*psiP(indxBdy(1):indxBdy(end)));
        f33=trapz(x(indxBdy(1):indxBdy(end)), psiH(indxBdy(1):indxBdy(end)));
        f53=trapz(x(indxBdy(1):indxBdy(end)),psiP(indxBdy(1):indxBdy(end)));
        f352=trapz(x(indxBdy(1):indxBdy(end)),(x(indxBdy(1):indxBdy(end))-Xc).*psiH(indxBdy(1):indxBdy(end)));
                          
        f11=f11d;
        f15=-(f151d+f152);
        f31=f31d;
        f35=-(f351d+f352); %% heave to roll
                  
        f51=-f51d;
        f55=f551d+f552;
      
        Ma11=real(f11); Ma15=real(f15);Ma13=real(f13); 
        Ma33=real(f33); Ma31=real(f31);Ma35=real(f35);
        Ma55=real(f55); Ma51=real(f51); Ma53=real(f53);
        Ma13=0;Ma31=0;Ma35=0;Ma53=0;
         
        
         % 
           
       
        B11=imag(f11);B13=imag(f13); B15=imag(f15);
        B31=imag(f31);B33=imag(f33); B35=imag(f35);
        B51=imag(f51);B53=imag(f53); B55=imag(f55);

        
%         dxi_Krad.z(1)     =  real(trapz(x(indl:indr)-Xc,(dxPsiH_x(indl:indr).^2+dzPsiH_x.^2)/2-dzPsiH_x.*nu.z(1,indl:indr)));%trapz(x(indl:indr)-Xc,(dxPsiH_x(indl:indr).^2-dzPsiH_x.^2)/2);%
%         
%         dxi_Krad.x(ii)= real(trapz(ZshipL,(dPhi1Sdx_zL.^2+dPhi1Sdz_zL.^2)./2-dPhi1Sdz_zL.*(dPhi1Sdx_zL.*draft./dx)./sqrt(1+(draft./dx)^2)) ...
%             -trapz(ZshipR,(dPhi3Sdx_zR.^2+dPhi3Sdz_zR.^2)./2+dPhi3Sdz_zR.*(dPhi3Sdx_zR.*draft./dx)./sqrt(1+(draft./dx)^2)));
%         
%         dxi_Krad.theta(1)  = real(-trapz(ZshipL,(ZshipL-Zc).*((dPhi1Pdx_zL.^2+dPhi1Pdz_zL.^2)./2-dPhi1Pdz_zL.*(dPhi1Pdx_zL.*draft./dx)./sqrt(1+(draft./dx)^2))) ...
%             +trapz(ZshipR,(ZshipL-Zc).*((dPhi3Pdx_zR.^2+dPhi3Pdz_zR.^2)./2+dPhi3Pdz_zR.*(dPhi3Pdx_zR.*draft./dx)./sqrt(1+(draft./dx)^2)))...
%             +trapz(x(indl:indr)-Xc,(x(indl:indr)-Xc-Xc).*((dxPsiP_x(indl:indr).^2+dzPsiP_x.^2)/2-dzPsiP_x.*nu.z(1,indl:indr))));
    end
    
    if strcmpi(Idmotion,'Surge')|| strcmpi(Idmotion,'Free')
        psiS=psiS.';%psiS(isnan(psiS))=0;
        FpsiEvS(indl)=FevLS;FpsiEvS(indr)=FevRS;
        
        sig1C=CoefS(Nm/2);
        FpsiS=(sig1C*dxtaug).*(D-max(dr));
        FpsiS=FpsiS.*chiSWL;
    end
    if strcmpi(Idmotion,'Heave')|| strcmpi(Idmotion,'Free')
        psiH=psiH.';%psiH(isnan(psiH))=0;
        FpsiEvH(indl)=FevLH;FpsiEvH(indr)=FevRH;
        sig1C=CoefH(Nm/2);
        FpsiH=(sig1C*dxtaug).*(D-max(dr));
        FpsiH=FpsiH.*chiSWL;
        
        FpsiH=FpsiH-(x-Xc)+dxD.*((D-draft)-((D-draft)/3-(x-Xc).^2./(D-draft))/2);
        FpsiH=FpsiH.*chiSWL;
    end
    if strcmpi(Idmotion,'Pitch')|| strcmpi(Idmotion,'Free')
        psiP=psiP.';%psiP(isnan(psiP))=0;
        FpsiEvP(indl)=FevLP;FpsiEvP(indr)=FevRP;
        sig1C=CoefP(Nm/2);
        FpsiP=(sig1C*dxtaug).*(D-max(dr));
        FpsiP=FpsiP.*chiSWL;
        
        FpsiP=FpsiP-(x-Xc).^2/2+(D-draft).^2/6+(dxD./(D-draft)).*((x-Xc).*((D-draft).^2/2)+(-(D-draft).^2/3.*(x-Xc)+(x-Xc).^3/3)./2);
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
    psi.surge=real(psiS);psi.heave=real(psiH);
    psi.pitch=real(psiP);
    Lpsiev.surge=LpsievS;Lpsiev.heave=LpsievH;
    Lpsiev.pitch=LpsievP;
    
    
