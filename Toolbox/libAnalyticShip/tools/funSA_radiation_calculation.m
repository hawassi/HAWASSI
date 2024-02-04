Idmotion='Free';
[Amat,BvectS,BvectH,BvectP] =fun_construct_matrices_radiation_Ndom(Ndom,...
                    0,0,tauGbdy,dxtauGbdy,Dbdy,dxDbdy,DZetabdy,dxDZetabdy,dr,...
                    Evmodes,kappaS,Idmotion,Xc,Zc,etaWl,dxetaWl,Xbdy);

             
                if strcmpi(Idmotion,'Surge')|| strcmpi(Idmotion,'Free')
                    CoefS=Amat\BvectS;
                    psiS=zeros(size(chiS));
                    Nm=Ndom*2-2;
                    idAd=0;idCoef1=3;idC=Ndom/2;idAdE=0;idCoef2=Nm-1;
                    eta1L=0;eta1R=0;
                    InQ0TEtaL=fun_IntegQn(kappaS(1,1),eta1L,Dbdy(1),-dr(1),eta1L);
                    InQ0TEtaR=fun_IntegQn(kappaS(1,2),eta1R,Dbdy(end),-dr(end),eta1R);
                    f11d=CoefS(1)*InQ0TEtaL-CoefS(Nm)*InQ0TEtaR;
                    
                    InQ0TEtaL=fun_IntegQn(kappaS(1,1),eta1L,Dbdy(1),-dr(1),eta1L);
                    InQ0TEtaR=fun_IntegQn(kappaS(1,2),eta1R,Dbdy(end),-dr(end),eta1R);
                    InZQ0TEtaL=fun_IntegZQn(kappaS(1,1),eta1L,Dbdy(1),-dr(1),eta1L);
                    InZQ0TEtaR=fun_IntegZQn(kappaS(1,2),eta1R,Dbdy(end),-dr(end),eta1R);
                    
                    f151d=-CoefS(1)*(InZQ0TEtaL-Zc*InQ0TEtaL)+CoefS(Nm)*(InZQ0TEtaR-Zc*InQ0TEtaR);
                    
                    for jj=1:(Ndom-2)/2
                        if jj==(Ndom-2)/2, idAdE=1; end
                        int1=indxBdy(jj)+idAd:indxBdy(jj+1);
                        psiS(int1)=psiS(int1)+CoefS(idCoef1).*tauGInterp(int1)+CoefS(idCoef1-1);
                        int2=indxBdy(end-jj)+idAdE:indxBdy(end-jj+1)-idAd;
                        psiS(int2)=psiS(int2)+CoefS(idCoef2).*tauGInterp(int2)+CoefS(idCoef2-1);
                        
                        if jj<(Ndom-2)/2
                            Intvarphi0Zetad_dprev=fun_IntegVarphin(0,Dbdy(jj+1),-dr(jj+1),-dr(jj));
                            f11d=f11d+(CoefS(idCoef1).*tauGInterp(indxBdy(jj+1))+CoefS(idCoef1-1))*Intvarphi0Zetad_dprev;
                            
                            Intvarphi0Zetad_dnext=fun_IntegVarphin(0,Dbdy(end-jj),-dr(end-jj),-dr(end-jj+1));
                            f11d=f11d-(CoefS(idCoef2).*tauGInterp(indxBdy(end-jj))+CoefS(idCoef2-1))*Intvarphi0Zetad_dnext;
                            
                            IntZvarphi0Zetad_dprev=fun_IntegZVarphin(0,Dbdy(jj+1),-dr(jj+1),-dr(jj));
                            IntZvarphi0Zetad_dnext=fun_IntegZVarphin(0,Dbdy(end-jj),-dr(end-jj),-dr(end-jj+1));
                            
                            f151d=f151d-(CoefS(idCoef1).*tauGInterp(indxBdy(jj+1))+CoefS(idCoef1-1))*(IntZvarphi0Zetad_dprev-Zc*Intvarphi0Zetad_dprev);
                            f151d=f151d+(CoefS(idCoef2).*tauGInterp(indxBdy(end-jj))+CoefS(idCoef2-1))*(IntZvarphi0Zetad_dnext-Zc*Intvarphi0Zetad_dnext);
                        end
                        idAd=1;idCoef1=idCoef1+2; idCoef2=idCoef2-2;
                    end
                    
                    psiS=psiS.*chiS;
                    
                    psiS(indxBdy(1))=CoefS(1);
                    psiS(indxBdy(end))=CoefS(Nm);
                    
                end
                
                
                if strcmpi(Idmotion,'Heave')|| strcmpi(Idmotion,'Free')
                    CoefH=Amat\BvectH;
                    GammaH(indx1:indx2)=((D+shipform(indx1:indx2)).^2-(xx(indx1:indx2)-Xc).^2)./2./(D+shipform(indx1:indx2));
                    % bb=GammaH(indx1:indx2);
                    psiH=GammaH.*chiSwl;
                    psiH(indx1)=GammaH(indx1+1);
                    psiH(indx2)=GammaH(indx2-1);
                    Nm=Ndom*2-2;
                    idAd=0;idCoef1=3;idC=Ndom/2;idAdE=0;idCoef2=Nm-1;
                    
                    eta1L=0;eta1R=0;
                    InQ0TEtaL=fun_IntegQn(kappaS(1,1),eta1L,Dbdy(1),-dr(1),eta1L);
                    InQ0TEtaR=fun_IntegQn(kappaS(1,2),eta1R,Dbdy(end),-dr(end),eta1R);
                    f31d=CoefH(1)*InQ0TEtaL-CoefH(Nm)*InQ0TEtaR;
                    
                    InZQ0TEtaL=fun_IntegZQn(kappaS(1,1),eta1L,Dbdy(1),-dr(1),eta1L);
                    InZQ0TEtaR=fun_IntegZQn(kappaS(1,2),eta1R,Dbdy(end),-dr(end),eta1R);
                    
                    f351d=-CoefH(1)*(InZQ0TEtaL-Zc*InQ0TEtaL)+CoefH(Nm)*(InZQ0TEtaR-Zc*InQ0TEtaR);
                    
                    for jj=1:(Ndom-2)/2
                        if jj==(Ndom-2)/2, idAdE=1; end
                        int1=indxBdy(jj)+idAd:indxBdy(jj+1);
                        psiH(int1)=psiH(int1)+CoefH(idCoef1).*tauGInterp(int1)+CoefH(idCoef1-1);
                        int2=indxBdy(end-jj)+idAdE:indxBdy(end-jj+1)-idAd;
                        psiH(int2)=psiH(int2)+CoefH(idCoef2).*tauGInterp(int2)+CoefH(idCoef2-1);
                        
                        
                        if jj<(Ndom-2)/2
                            Intvarphi0Zetad_dprev=fun_IntegVarphin(0,Dbdy(jj+1),-dr(jj+1),-dr(jj));
                            f31d=f31d+(CoefH(idCoef1).*tauGInterp(indxBdy(jj+1))+CoefH(idCoef1-1))*Intvarphi0Zetad_dprev;
                            
                            Intvarphi0Zetad_dnext=fun_IntegVarphin(0,Dbdy(end-jj),-dr(end-jj),-dr(end-jj+1));
                            f31d=f31d-(CoefH(idCoef2).*tauGInterp(indxBdy(end-jj))+CoefH(idCoef2-1))*Intvarphi0Zetad_dnext;
                            
                            IntZvarphi0Zetad_dprev=fun_IntegZVarphin(0,Dbdy(jj+1),-dr(jj+1),-dr(jj));
                            IntZvarphi0Zetad_dnext=fun_IntegZVarphin(0,Dbdy(end-jj),-dr(end-jj),-dr(end-jj+1));
                            
                            f351d=f351d-(CoefH(idCoef1).*tauGInterp(indxBdy(jj+1))+CoefH(idCoef1-1))*(IntZvarphi0Zetad_dprev-Zc*Intvarphi0Zetad_dprev);
                            f351d=f351d+(CoefH(idCoef2).*tauGInterp(indxBdy(end-jj))+CoefH(idCoef2-1))*(IntZvarphi0Zetad_dnext-Zc*Intvarphi0Zetad_dnext);
                        end
                        
                        idAd=1;idCoef1=idCoef1+2; idCoef2=idCoef2-2;
                    end
                    psiH=psiH.*chiSwl;
                end
                
                if strcmpi(Idmotion,'Pitch')|| strcmpi(Idmotion,'Free')
                    CoefP=Amat\BvectP;
                    
                    GammaP(indx1:indx2)=((D+shipform(indx1:indx2)).^2.*(xx(indx1:indx2)-Xc)-(xx(indx1:indx2)-Xc).^3/3)./2./(D+shipform(indx1:indx2));
                    
                    psiP=GammaP.*chiSwl;
                    psiP(indx1)=GammaP(indx1+1);
                    psiP(indx2)=GammaP(indx2-1);
                    %
                    Nm=Ndom*2-2;
                    idAd=0;idCoef1=3;idC=Ndom/2;idAdE=0;idCoef2=Nm-1;
                    eta1L=0;eta1R=0;
                    InQ0TEtaL=fun_IntegQn(kappaS(1,1),eta1L,Dbdy(1),-dr(1),eta1L);
                    InQ0TEtaR=fun_IntegQn(kappaS(1,2),eta1R,Dbdy(end),-dr(end),eta1R);
                    InZQ0TEtaL=fun_IntegZQn(kappaS(1,1),eta1L,Dbdy(1),-dr(1),eta1L);
                    InZQ0TEtaR=fun_IntegZQn(kappaS(1,2),eta1R,Dbdy(end),-dr(end),eta1R);
                    f51d=CoefP(1)*InQ0TEtaL-CoefP(Nm)*InQ0TEtaR;
                    
                    f551d=-CoefP(1)*(InZQ0TEtaL-Zc*InQ0TEtaL)+CoefP(Nm)*(InZQ0TEtaR-Zc*InQ0TEtaR);
                    
                    for jj=1:(Ndom-2)/2
                        if jj==(Ndom-2)/2, idAdE=1; end
                        int1=indxBdy(jj)+idAd:indxBdy(jj+1);
                        psiP(int1)=psiP(int1)+CoefP(idCoef1).*tauGInterp(int1)+CoefP(idCoef1-1);
                        int2=indxBdy(end-jj)+idAdE:indxBdy(end-jj+1)-idAd;
                        psiP(int2)=psiP(int2)+CoefP(idCoef2).*tauGInterp(int2)+CoefP(idCoef2-1);
                        if jj<(Ndom-2)/2
                            Intvarphi0Zetad_dprev=fun_IntegVarphin(0,Dbdy(jj+1),-dr(jj+1),-dr(jj));
                            IntZvarphi0Zetad_dprev=fun_IntegZVarphin(0,Dbdy(jj+1),-dr(jj+1),-dr(jj));
                            IntGammaPZZ0_dprev=fun_IntGammaPZZ0(Dbdy(jj+1),Xbdy(jj+1)-Xc,Zc,-dr(jj+1),-dr(jj));
                            f551d=f551d-IntGammaPZZ0_dprev-(CoefP(idCoef1).*tauGInterp(indxBdy(jj+1))+CoefP(idCoef1-1))*(IntZvarphi0Zetad_dprev-Zc*Intvarphi0Zetad_dprev);
                            f51d=f51d+(CoefP(idCoef1).*tauGInterp(indxBdy(jj+1))+CoefP(idCoef1-1))*Intvarphi0Zetad_dprev;
                            
                            Intvarphi0Zetad_dnext=fun_IntegVarphin(0,Dbdy(end-jj),-dr(end-jj),-dr(end-jj+1));
                            IntZvarphi0Zetad_dnext=fun_IntegZVarphin(0,Dbdy(end-jj),-dr(end-jj),-dr(end-jj+1));
                            IntGammaPZZ0_dnext=fun_IntGammaPZZ0(Dbdy(end-jj),Xbdy(end-jj)-Xc,Zc,-dr(end-jj),-dr(end-jj+1));
                            f551d=f551d+IntGammaPZZ0_dnext+(CoefP(idCoef2).*tauGInterp(indxBdy(end-jj))+CoefP(idCoef2-1))*(IntZvarphi0Zetad_dnext-Zc*Intvarphi0Zetad_dnext);
                            f51d=f51d-(CoefP(idCoef2).*tauGInterp(indxBdy(end-jj))+CoefP(idCoef2-1))*Intvarphi0Zetad_dnext;
                        end
                        
                        idAd=1;idCoef1=idCoef1+2; idCoef2=idCoef2-2;
                    end
                    psiP=psiP.*chiSwl;
                end
                
                if Evmodes>0
                    for m=1:Evmodes
                        mu_m=m*pi./(DZeta);
                        
                        InQmTEtaL=fun_IntegQn(kappaS(m+1,1),eta1L,Dbdy(1),-dr(1),eta1L);
                        InQmTEtaR=fun_IntegQn(kappaS(m+1,2),eta1R,Dbdy(end),-dr(end),eta1R);
                        
                        
                        if strcmpi(Idmotion,'Surge')|| strcmpi(Idmotion,'Free')
                            psiSev=zeros(size(xx));
                            %
                            idAd=0;idAdE=0; idCoef1=Nm*m+2; idCoef2=Nm*m+Nm-2;
                            f11d=f11d+CoefS(Nm*m+1)*InQmTEtaL-CoefS(Nm*(m+1))*InQmTEtaR;
                            
                            InZQmTEtaL=fun_IntegZQn(kappaS(m+1,1),eta1L,Dbdy(1),-dr(1),eta1L);
                            InZQmTEtaR=fun_IntegZQn(kappaS(m+1,2),eta1R,Dbdy(end),-dr(end),eta1R);
                            f151d=f151d-CoefS(Nm*m+1)*(InZQmTEtaL-Zc*InQmTEtaL)+CoefS(Nm*(m+1))*(InZQmTEtaR-Zc*InQmTEtaR);
                            
                            for jj=1:(Ndom-2)/2
                                if jj==(Ndom-2)/2, idAdE=1; end
                                int1=indxBdy(jj)+idAd:indxBdy(jj+1);
                                psiSev(int1)=((CoefS(idCoef1).*exp(mu_m(int1).*(xx(int1)-Xbdy(jj)))...
                                    +CoefS(idCoef1+1).*exp(-mu_m(int1).*(xx(int1)-Xbdy(jj+1)))).*cos(mu_m(int1).*(DZeta(int1))));
                                int2=indxBdy(end-jj)+idAdE:indxBdy(end-jj+1)-idAd;
                                psiSev(int2)=((CoefS(idCoef2).*exp(mu_m(int2).*(xx(int2)-Xbdy(end-jj)))...
                                    +CoefS(idCoef2+1).*exp(-mu_m(int2).*(xx(int2)-Xbdy(end-jj+1)))).*cos(mu_m(int2).*(DZeta(int2))));
                                
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
                                end
                                
                                idAd=1;idCoef1=idCoef1+2;idCoef2=idCoef2-2;
                            end
                            psiSev=psiSev.*chiS;
                            psiSev(indxBdy(1))=CoefS(Nm*m+1);
                            psiSev(indxBdy(end))=CoefS(Nm*(m+1));
                            psiS=psiS+psiSev;
                        end
                        
                        if strcmpi(Idmotion,'Heave')|| strcmpi(Idmotion,'Free')
                            psiHev=zeros(size(xx));
                            %
                            idAd=0;idAdE=0; idCoef1=Nm*m+2; idCoef2=Nm*m+Nm-2;
                            
                            
                            f31d=f31d+CoefH(Nm*m+1)*InQmTEtaL-CoefH(Nm*(m+1))*InQmTEtaR;
                            
                            InZQmTEtaL=fun_IntegZQn(kappaS(m+1,1),eta1L,Dbdy(1),-dr(1),eta1L);
                            InZQmTEtaR=fun_IntegZQn(kappaS(m+1,2),eta1R,Dbdy(end),-dr(end),eta1R);
                            f351d=f351d-CoefH(Nm*m+1)*(InZQmTEtaL-Zc*InQmTEtaL)+CoefH(Nm*(m+1))*(InZQmTEtaR-Zc*InQmTEtaR);
                            
                            
                            for jj=1:(Ndom-2)/2
                                if jj==(Ndom-2)/2, idAdE=1; end
                                int1=indxBdy(jj)+idAd:indxBdy(jj+1);
                                psiHev(int1)=((CoefH(idCoef1).*exp(mu_m(int1).*(xx(int1)-Xbdy(jj)))...
                                    +CoefH(idCoef1+1).*exp(-mu_m(int1).*(xx(int1)-Xbdy(jj+1)))).*cos(mu_m(int1).*(DZeta(int1))));
                                int2=indxBdy(end-jj)+idAdE:indxBdy(end-jj+1)-idAd;
                                psiHev(int2)=((CoefH(idCoef2).*exp(mu_m(int2).*(xx(int2)-Xbdy(end-jj)))...
                                    +CoefH(idCoef2+1).*exp(-mu_m(int2).*(xx(int2)-Xbdy(end-jj+1)))).*cos(mu_m(int2).*(DZeta(int2))));
                                
                                
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
                            psiHev=psiHev.*chiSwl;
                            psiH=psiH+psiHev;
                        end
                        
                        if strcmpi(Idmotion,'Pitch')|| strcmpi(Idmotion,'Free')
                            psiPev=zeros(size(xx));
                            
                            InQmTEtaL=fun_IntegQn(kappaS(m+1,1),eta1L,Dbdy(1),-dr(1),eta1L);
                            InQmTEtaR=fun_IntegQn(kappaS(m+1,2),eta1R,Dbdy(end),-dr(end),eta1R);
                            InZQmTEtaL=fun_IntegZQn(kappaS(m+1,1),eta1L,Dbdy(1),-dr(1),eta1L);
                            InZQmTEtaR=fun_IntegZQn(kappaS(m+1,2),eta1R,Dbdy(end),-dr(end),eta1R);
                            
                            f51d=f51d+CoefP(Nm*m+1)*InQmTEtaL-CoefP(Nm*(m+1))*InQmTEtaR;
                            
                            f551d=f551d-CoefP(Nm*m+1)*(InZQmTEtaL-Zc*InQmTEtaL)+CoefP(Nm*(m+1))*(InZQmTEtaR-Zc*InQmTEtaR);
                            
                            idAd=0;idAdE=0; idCoef1=Nm*m+2; idCoef2=Nm*m+Nm-2;
                            for jj=1:(Ndom-2)/2
                                if jj==(Ndom-2)/2, idAdE=1; end
                                int1=indxBdy(jj)+idAd:indxBdy(jj+1);
                                psiPev(int1)=((CoefP(idCoef1).*exp(mu_m(int1).*(xx(int1)-Xbdy(jj)))...
                                    +CoefP(idCoef1+1).*exp(-mu_m(int1).*(xx(int1)-Xbdy(jj+1)))).*cos(mu_m(int1).*(DZeta(int1))));
                                int2=indxBdy(end-jj)+idAdE:indxBdy(end-jj+1)-idAd;
                                psiPev(int2)=((CoefP(idCoef2).*exp(mu_m(int2).*(xx(int2)-Xbdy(end-jj)))...
                                    +CoefP(idCoef2+1).*exp(-mu_m(int2).*(xx(int2)-Xbdy(end-jj+1)))).*cos(mu_m(int2).*(DZeta(int2))));
                                
                                
                                if jj<(Ndom-2)/2
                                    mu_mLd=m*pi./(Dbdy(jj)-dr(jj));
                                    mu_mRd=m*pi./(Dbdy(end-jj+1)-dr(end-jj+1));
                                    IntvarphimZetad_dprev=fun_IntegVarphin(mu_mLd,Dbdy(jj+1),-dr(jj+1),-dr(jj));
                                    IntZvarphimZetad_dprev=fun_IntegZVarphin(mu_mLd,Dbdy(jj+1),-dr(jj+1),-dr(jj));
                                    IntvarphimZetad_dnext=fun_IntegVarphin(mu_mRd,Dbdy(end-jj),-dr(end-jj),-dr(end-jj+1));
                                    IntZvarphimZetad_dnext=fun_IntegZVarphin(mu_mRd,Dbdy(end-jj),-dr(end-jj),-dr(end-jj+1));
                                    
                                    f51d=f51d+(CoefP(idCoef1).*exp(mu_mLd.*(Xbdy(jj+1)-Xbdy(jj)))...
                                        +CoefP(idCoef1+1)).*IntvarphimZetad_dprev;
                                    f51d=f51d-(CoefP(idCoef2)+CoefP(idCoef2+1).*exp(-mu_mRd.*(Xbdy(end-jj)-Xbdy(end-jj+1))))*IntvarphimZetad_dnext;
                                    
                                    f551d=f551d-(CoefP(idCoef1).*exp(mu_mLd.*(Xbdy(jj+1)-Xbdy(jj)))...
                                        +CoefP(idCoef1+1)).*(IntZvarphimZetad_dprev-Zc*IntvarphimZetad_dprev);
                                    f551d=f551d+(CoefP(idCoef2)+CoefP(idCoef2+1).*...
                                        exp(-mu_mRd.*(Xbdy(end-jj)-Xbdy(end-jj+1))))*(IntZvarphimZetad_dnext-Zc*IntvarphimZetad_dnext);
                                end
                                idAd=1;idCoef1=idCoef1+2;idCoef2=idCoef2-2;
                            end
                            psiPev=psiPev.*chiSwl;
                            psiP=psiP+psiPev;
                        end
                        
                    end
                end
                
                    psiSXL=interp1(xx(indxBdy(1):indxBdy(end)),psiS(indxBdy(1):indxBdy(end)),XshipL);
                    psiSXR=interp1(xx(indxBdy(1):indxBdy(end)),psiS(indxBdy(1):indxBdy(end)),XshipR);
                    %                 f11x=-trapz(XshipL,-psiSXL.*dxZshipL)+trapz(XshipR,-psiSXR.*dxZshipR);
                    %                 f11=trapz(ZshipL,psiSXL)-trapz(ZshipR,psiSXR);
                    f11=f11d;
                    f152=trapz(xx(indxBdy(1):indxBdy(end)),(xx(indxBdy(1):indxBdy(end))-Xc).*psiS(indxBdy(1):indxBdy(end)));
                    f15=-(f151d+f152); %% sway to roll
                    f13= trapz(xx(indxBdy(1):indxBdy(end)),psiS(indxBdy(1):indxBdy(end)));
                    
                    psiS=real(psiS);
                  
                
                    f31=f31d;
                    f33=trapz(xx(indxBdy(1):indxBdy(end)),psiH(indxBdy(1):indxBdy(end)));
                    f352=trapz(xx(indxBdy(1):indxBdy(end)),(xx(indxBdy(1):indxBdy(end))-Xc).*psiH(indxBdy(1):indxBdy(end)));
                    f35=-(f351d+f352); %% heave to roll
                    
                    psiH=real(psiH);
                    %      if Idplot==1
                    %      figure;
                  
                
                     psiPXL=interp1(xx(indxBdy(1):indxBdy(end)),psiP(indxBdy(1):indxBdy(end)),XshipL);
                    psiPXR=interp1(xx(indxBdy(1):indxBdy(end)),psiP(indxBdy(1):indxBdy(end)),XshipR);
                    f51=-f51d;
                    f53=trapz(xx(indxBdy(1):indxBdy(end)),psiP(indxBdy(1):indxBdy(end)));
                    f55x=-trapz(XshipL,psiPXL.*((ZshipL-Zc).*dxZshipL+(XshipL-Xc)))+trapz(XshipR,psiPXR.*((ZshipR-Zc).*dxZshipR+(XshipR-Xc)));
                    f55z=-trapz(ZshipL,psiPXL.*((ZshipL-Zc)+dzX_L.*(XshipL-Xc)))+trapz(ZshipR,psiPXR.*((ZshipR-Zc)+dzX_R.*(XshipR-Xc)));
                    f552=trapz(xx(indxBdy(1):indxBdy(end)),(xx(indxBdy(1):indxBdy(end))-Xc).*psiP(indxBdy(1):indxBdy(end)));
                    f55=f551d+f552;
                    
                    
                    psiP=real(psiP);
                    %      if Idplot==1
                    %      figure;
                    
%                     subplot(3,1,1)
%                     plot(xx,psiS,'r');
%                     title(['Added Mass: ',num2str(real(f11))])
%                     subplot(3,1,2)
%                     plot(xx,psiH,'r');
%                     title(['Added Mass: ',num2str(real(f33))])
%                     subplot(3,1,3)
%                     plot(xx,psiP,'r');
%                     title(['Added Mass: ',num2str(real(f55x))])
%                     pause(0.01);
              
                rad_waveAmpl(ii,1)=w0(ii);
                rad_waveAmpl(ii,2:3)=abs(-((1i.*w0(ii))^2.*[CoefS(1) CoefS(Nm)]/g));
                rad_waveAmpl(ii,4:5)=abs(-((1i.*w0(ii))^2.*[CoefH(1) CoefH(Nm)]/g));
                rad_waveAmpl(ii,6:7)=abs(-((1i.*w0(ii))^2.*[CoefP(1) CoefP(Nm)]/g));
                rad_addedMass(ii,1)=w0(ii);
                rad_addedMass(ii,2:4) =[real(f11)  real(f13) real(f15)];
                rad_addedMass(ii,5:7) =[real(f31)  real(f33) real(f35)];
                rad_addedMass(ii,8:10)=[real(f51) real(f53) real(f55)];
                rad_dampCoef(ii,1)   =w0(ii);
                rad_dampCoef(ii,2:4) =w0(ii).*[imag(f11) imag(f13) imag(f15)];
                rad_dampCoef(ii,5:7) =w0(ii).*[imag(f31) imag(f33) imag(f35)];
                rad_dampCoef(ii,8:10)=w0(ii).*[imag(f51) imag(f53) imag(f55)];
                    
                  
               