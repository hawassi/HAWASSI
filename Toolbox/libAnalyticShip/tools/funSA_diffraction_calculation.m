 [AmatD,BvectD] =funSA_construct_matrices_diffraction_Ndom(Ndom,...
                    0,0,tauGbdy,dxtauGbdy,Dbdy,dxDbdy,DZetabdy,dxDZetabdy,dr,...
                    Evmodes,kappaS,etaWl,dxetaWl,Xbdy,IncAmpl);
                CoefD=AmatD\BvectD;
             
                
                phiF=zeros(size(chiS));
                Nm=Ndom*2-2;
                
                
                idAd=0;idCoef1=3;idC=Ndom/2;idAdE=0;idCoef2=Nm-1;
                eta1L=0;eta1R=0;
                InQ0TEtaL=fun_IntegQn(kappaS(1,1),eta1L,Dbdy(1),-dr(1),eta1L);
                InQ0TEtaR=fun_IntegQn(kappaS(1,2),eta1R,Dbdy(end),-dr(end),eta1R);
                betaxd=(IncAmpl+CoefD(1))*InQ0TEtaL-CoefD(Nm)*InQ0TEtaR;
                
                InQ0TEtaL=fun_IntegQn(kappaS(1,1),eta1L,Dbdy(1),-dr(1),eta1L);
                InQ0TEtaR=fun_IntegQn(kappaS(1,2),eta1R,Dbdy(end),-dr(end),eta1R);
                InZQ0TEtaL=fun_IntegZQn(kappaS(1,1),eta1L,Dbdy(1),-dr(1),eta1L);
                InZQ0TEtaR=fun_IntegZQn(kappaS(1,2),eta1R,Dbdy(end),-dr(end),eta1R);
                
                betathetad=-(IncAmpl+CoefD(1))*(InZQ0TEtaL-Zc*InQ0TEtaL)+CoefD(Nm)*(InZQ0TEtaR-Zc*InQ0TEtaR);
                
                
                
                for jj=1:(Ndom-2)/2
                    if jj==(Ndom-2)/2, idAdE=1; end
                    int1=indxBdy(jj)+idAd:indxBdy(jj+1);
                    phiF(int1)=phiF(int1)+CoefD(idCoef1).*tauGInterp(int1)+CoefD(idCoef1-1);
                    int2=indxBdy(end-jj)+idAdE:indxBdy(end-jj+1)-idAd;
                    phiF(int2)=phiF(int2)+CoefD(idCoef2).*tauGInterp(int2)+CoefD(idCoef2-1);
                    
                    if jj<(Ndom-2)/2
                        Intvarphi0Zetad_dprev=fun_IntegVarphin(0,Dbdy(jj+1),-dr(jj+1),-dr(jj));
                        betaxd=betaxd+(CoefD(idCoef1).*tauGInterp(indxBdy(jj+1))+CoefD(idCoef1-1))*Intvarphi0Zetad_dprev;
                        
                        Intvarphi0Zetad_dnext=fun_IntegVarphin(0,Dbdy(end-jj),-dr(end-jj),-dr(end-jj+1));
                        betaxd=betaxd-(CoefD(idCoef2).*tauGInterp(indxBdy(end-jj))+CoefD(idCoef2-1))*Intvarphi0Zetad_dnext;
                        
                        IntZvarphi0Zetad_dprev=fun_IntegZVarphin(0,Dbdy(jj+1),-dr(jj+1),-dr(jj));
                        betathetad=betathetad-(CoefD(idCoef1).*tauGInterp(indxBdy(jj+1))+CoefD(idCoef1-1))*(IntZvarphi0Zetad_dprev-Zc*Intvarphi0Zetad_dprev);
                        
                        IntZvarphi0Zetad_dnext=fun_IntegZVarphin(0,Dbdy(end-jj),-dr(end-jj),-dr(end-jj+1));
                        betathetad=betathetad+(CoefD(idCoef2).*tauGInterp(indxBdy(end-jj))+CoefD(idCoef2-1))*(IntZvarphi0Zetad_dnext-Zc*Intvarphi0Zetad_dnext);
                    end
                    idAd=1;idCoef1=idCoef1+2; idCoef2=idCoef2-2;
                end
                phiF=phiF.*chiS;
                
                %                 phiF(indxBdy(1))=IncAmpl+CoefD(1);
                %                 phiF(indxBdy(end))=CoefD(Nm);
                %
                if Evmodes>0
                    for m=1:Evmodes
                        mu_m=m*pi./(DZeta);
                        phiFev=zeros(size(xx));
                        %
                        idAd=0;idAdE=0; idCoef1=Nm*m+2; idCoef2=Nm*m+Nm-2;
                        InQmTEtaL=fun_IntegQn(kappaS(m+1,1),eta1L,Dbdy(1),-dr(1),eta1L);
                        InQmTEtaR=fun_IntegQn(kappaS(m+1,2),eta1R,Dbdy(end),-dr(end),eta1R);
                        betaxd=betaxd+CoefD(Nm*m+1)*InQmTEtaL-CoefD(Nm*(m+1))*InQmTEtaR;
                        InZQmTEtaL=fun_IntegZQn(kappaS(m+1,1),eta1L,Dbdy(1),-dr(1),eta1L);
                        InZQmTEtaR=fun_IntegZQn(kappaS(m+1,2),eta1R,Dbdy(end),-dr(end),eta1R);
                        
                        betathetad=betathetad-CoefD(Nm*m+1)*(InZQmTEtaL-Zc*InQmTEtaL)+CoefD(Nm*(m+1))*(InZQmTEtaR-Zc*InQmTEtaR);
                        
                        for jj=1:(Ndom-2)/2
                            if jj==(Ndom-2)/2, idAdE=1; end
                            int1=indxBdy(jj)+idAd:indxBdy(jj+1);
                            phiFev(int1)=((CoefD(idCoef1).*exp(mu_m(int1).*(xx(int1)-Xbdy(jj)))...
                                +CoefD(idCoef1+1).*exp(-mu_m(int1).*(xx(int1)-Xbdy(jj+1)))).*cos(mu_m(int1).*(DZeta(int1))));
                            int2=indxBdy(end-jj)+idAdE:indxBdy(end-jj+1)-idAd;
                            phiFev(int2)=((CoefD(idCoef2).*exp(mu_m(int2).*(xx(int2)-Xbdy(end-jj)))...
                                +CoefD(idCoef2+1).*exp(-mu_m(int2).*(xx(int2)-Xbdy(end-jj+1)))).*cos(mu_m(int2).*(DZeta(int2))));
                            
                            if jj<(Ndom-2)/2
                                mu_mLd=m*pi./(D-dr(jj));
                                mu_mRd=m*pi./(D-dr(end-jj+1));
                                IntvarphimZetad_dprev=fun_IntegVarphin(mu_mLd,Dbdy(jj+1),-dr(jj+1),-dr(jj));
                                betaxd=betaxd+(CoefD(idCoef1).*exp(mu_mLd.*(Xbdy(jj+1)-Xbdy(jj)))...
                                    +CoefD(idCoef1+1)).*IntvarphimZetad_dprev;
                                IntvarphimZetad_dnext=fun_IntegVarphin(mu_mRd,Dbdy(end-jj),-dr(end-jj),-dr(end-jj+1));
                                betaxd=betaxd-(CoefD(idCoef2)+CoefD(idCoef2+1).*exp(-mu_mRd.*(Xbdy(end-jj)-Xbdy(end-jj+1))))*IntvarphimZetad_dnext;
                                
                                IntZvarphimZetad_dprev=fun_IntegZVarphin(mu_mLd,Dbdy(jj+1),-dr(jj+1),-dr(jj));
                                betathetad=betathetad-(CoefD(idCoef1).*exp(mu_mLd.*(Xbdy(jj+1)-Xbdy(jj)))...
                                    +CoefD(idCoef1+1)).*(IntZvarphimZetad_dprev-Zc*IntvarphimZetad_dprev);
                                
                                IntZvarphimZetad_dnext=fun_IntegZVarphin(mu_mRd,Dbdy(end-jj),-dr(end-jj),-dr(end-jj+1));
                                
                                betathetad=betathetad+(CoefD(idCoef2)+CoefD(idCoef2+1).*...
                                    exp(-mu_mRd.*(Xbdy(end-jj)-Xbdy(end-jj+1))))*(IntZvarphimZetad_dnext-Zc*IntvarphimZetad_dnext);
                                
                            end
                            
                            idAd=1;idCoef1=idCoef1+2;idCoef2=idCoef2-2;
                        end
                        phiFev=phiFev.*chiS;
                        phiFev(indxBdy(1))=CoefD(Nm*m+1);
                        phiFev(indxBdy(end))=CoefD(Nm*(m+1));
                        phiF=phiF+phiFev;
                        
                    end
                end
                
                phiFXL=interp1(xx(indxBdy(1):indxBdy(end)),phiF(indxBdy(1):indxBdy(end)),XshipL);
                phiFXR=interp1(xx(indxBdy(1):indxBdy(end)),phiF(indxBdy(1):indxBdy(end)),XshipR);
                betaxx=-trapz(XshipL,-phiFXL.*dxZshipL)+trapz(XshipR,-phiFXR.*dxZshipR);
                betaxz=trapz(ZshipL,phiFXL)-trapz(ZshipR,phiFXR);
                betazm=trapz(xx(indxBdy(1):indxBdy(end)),phiF(indxBdy(1):indxBdy(end)));
                betathetad=betathetad+...
                    trapz(xx(indxBdy(1):indxBdy(end)),(xx(indxBdy(1):indxBdy(end))-Xc).*phiF(indxBdy(1):indxBdy(end)));
                
                
                betax=-1i*g./w0(ii)*betaxd; % we change amplitude of phi according amplitude of eta
                betaz=-1i*g./w0(ii)*betazm;
                betatheta=-1i*g./w0(ii)*(betathetad);
                
%                 if IDTestcase==13
%                     factx=g*IncAmpl*draft;facttheta=g*IncAmpl*B*draft;factz=g*IncAmpl*B;
%                 elseif IDTestcase==14 || IDTestcase==15 || IDTestcase==16
%                     factx=g*IncAmpl*B;facttheta=g*IncAmpl*B*b;factz=g*IncAmpl*B;
%                 elseif IDTestcase==18 || IDTestcase==19
%                     factx=g*b*IncAmpl;factz=g*b*IncAmpl;facttheta=g*IncAmpl*b^2;
%                 elseif IDTestcase==20
%                     factx=g*D*IncAmpl;factz=g*D*IncAmpl;facttheta=g*IncAmpl*D^2;
%                 else
%                    factx=g*AreaO*kAinc;factz=g*B*IncAmpl;facttheta=g*B^3*kAinc/12;
%                 end
                absFx=abs((-1i.*w0(ii).*betax));
                absFz=abs((-1i.*w0(ii).*betaz));
                absFtheta=abs((-1i.*w0(ii).*betatheta));
               
                
                phiF=real(phiF);
                dtt=0.1;
                t=[0:dtt:600];
                wwt=freqspace(t.');
%                 subplot(2,1,1)
%                 plot(xx,phiF,'r',xx,shipform_aprox,'k');
%                 title(['beta x: ',num2str(real(betaxx))])
               
                Fxamp=-1i.*w0(ii).*betax;
                Fzamp=-1i.*w0(ii).*betaz;
                Fthetaamp=-1i.*w0(ii).*betatheta;
                eta_t=real(IncAmpl.*exp(-1i.*w0(ii).*t));
                phi_t=real(-1i*g./w0(ii)*IncAmpl.*exp(-1i.*w0(ii).*t));
                Fx_t=real(Fxamp.*exp(-1i.*w0(ii).*t));%gradient(real(betax(ii).*exp(-1i.*w0(ii).*t)),dtt);
                Fz_t=real(Fzamp.*exp(-1i.*w0(ii).*t));%gradient(real(betaz(ii).*exp(-1i.*w0(ii).*t)),dtt);
                Ftheta_t=real(Fthetaamp.*exp(-1i.*w0(ii).*t));%gradient(real(betatheta(ii).*exp(-1i.*w0(ii).*t)),dtt);
                     
                phasePhi2eta=acos(fun_corr(eta_t.',phi_t.'));
                phaseFx2eta=acos(fun_corr(eta_t.',Fx_t.'));
                phaseFz2eta=acos(fun_corr(eta_t.',Fz_t.'));
                phaseFtheta2eta=acos(fun_corr(eta_t.',Ftheta_t.'));
                
%                 subplot(2,1,2)
%                 plot(t,eta_t./max(eta_t),'b',t,phi_t./max(phi_t),'--b',t,Fx_t./max(Fx_t),'--r',t,Fz_t./max(Fz_t),'--g',t,Ftheta_t./max(Ftheta_t),'-.c');
%                 legend('eta','phi','Fx','Fz','Ftheta')
%                 pause(0.01)%
               
                Rcoef=abs(CoefD(1)/IncAmpl);Tcoef=abs(CoefD(Nm)/IncAmpl);
               
                diff_RefTransCoef(ii,:)=[w0(ii) Rcoef Tcoef];
                diff_Forces(ii,:)=[w0(ii) absFx absFz absFtheta phaseFx2eta phaseFz2eta phaseFtheta2eta];