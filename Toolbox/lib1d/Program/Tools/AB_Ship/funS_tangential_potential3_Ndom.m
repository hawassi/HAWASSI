function [phif,phifZ0,Fs,Fev,Lev,beta,dxi_K,etaWl,dxetaWl,tauAz]=funS_tangential_potential3_Ndom(calclin,k,x,xShip,phi,eta,chiship,chiwaterline,depth,Evmodes,XZrot0,kappaS,taug,shipformOnGrid,shapedat)
Nship=1;
%%%% find anzatz tangential potential (before apply impermeability condition)
chiSWl=(chiship(:,end)+chiwaterline(:,end));
%tauAz=fun_tangential_potential_anzats(k,x,L0,F0,taug,eta,phi,Nonlin,ChiAdj,xShip,chiSWl,Nship);
tauAz=zeros(size(chiSWl));%
%%%%


Fs=zeros(size(phi));
phif=phi;

beta.z       =zeros(Nship,1);
beta.x       =zeros(Nship,1);
beta.theta   =zeros(Nship,1);

dxphi=ifft(1i.*k.*fft(phi),'symmetric');
dxeta=ifft(1i.*k.*fft(eta),'symmetric');
dxi_K.z       =zeros(Nship,1);
dxi_K.x       =zeros(Nship,1);
dxi_K.theta   =zeros(Nship,1);
dx=x(2)-x(1);
D=depth;
dxD=gradient(D,dx);

Xwl=[xShip(1,1) xShip(1,3)];
indWl=closest(x,Xwl);
Dwl=D(indWl);
dxDwl=dxD(indWl);

Zc=XZrot0(1,2);%%
Xc=xShip(1,2);%
kappaS(1,:)=-kappaS(1,:);

indl=indWl(1);indr=indWl(2);

    phiL=(phi(indl)+phi(indl-1)+phi(indl-2)+phi(indl-3))/4;%
    phiR=(phi(indr+3)+phi(indr+2)+phi(indr+1)+phi(indr))/4;%
    dxphiL=(dxphi(indl)+dxphi(indl-1)+dxphi(indl-2)+dxphi(indl-3))/4;%
    dxphiR=(dxphi(indr+3)+dxphi(indr+2)+dxphi(indr+1)+dxphi(indr))/4;%
    etaL=0;%(eta(indl)+eta(indl-1)+eta(indl-2)+eta(indl-3))/4;%
    etaR=0;%(eta(indr+3)+eta(indr+2)+eta(indr+1)+eta(indr))/4;%
    dxetaL=0;%(dxeta(indl)+dxeta(indl-1)+dxeta(indl-2)+dxeta(indl-3))/4;%
    dxetaR=0;%(dxeta(indr+3)+dxeta(indr+2)+dxeta(indr+1)+dxeta(indr))/4;%   
%        if calclin==1
%             etaL=0;etaR=0;
%       end
        phiWl=[phiL phiR];
        etaWl=[etaL etaR];
      
        dxphiWl=[dxphiL dxphiR];
        dxetaWl=[dxetaL dxetaR];
        
        tauAzWl=tauAz(indWl);
        taugWl=taug(indWl);
        dxtauAz=gradient(tauAz,dx);
        dxtaug=gradient(taug,dx);
%         dxtaug(indWl(1))=dxtaug(indWl(1)+1);
%         dxtaug(indWl(2))=dxtaug(indWl(2)-1);
%         dxtauAz(indWl(1))=dxtauAz(indWl(1)+1);
%         dxtauAz(indWl(2))=dxtauAz(indWl(2)-1);
        
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         dxtauAzWl=dxtauAz(indWl);
%         dxtaugWl=dxtaug(indWl);
        Ndom=shapedat(1).Ndom;
%         XshipL=shapedat(1).XshipL;
%         XshipR=shapedat(1).XshipR;
        ZshipL=shapedat(1).ZshipL;
        ZshipR=shapedat(1).ZshipR;
%         dxZshipL=shapedat(1).dxZshipL;
%         dxZshipR=shapedat(1).dxZshipR;
        
        indxBdy=shapedat(1).indxBdy;
        
        Xbdy=shapedat(1).Xbdy;
        dr=shapedat(1).dr_approx;
       
        Dbdy=D(indxBdy);
        dxDbdy=dxD(indxBdy);
        tauGbdy=taug(indxBdy);
        dxtauGbdy=dxtaug(indxBdy);
        
        
        DZeta=D+shipformOnGrid(end,:).';
        %dxDZeta=funOprt_FDGradient1d(DZeta,dx,2);
        DZetabdy=DZeta(indxBdy);%[DZeta(indx1) DZeta(indxlc) DZeta(indc) DZeta(indxcr) DZeta(indx2)];
        dxDZetabdy=zeros(size(DZetabdy));%%dxDZeta(indxBdy);%[dxDbdy(1) dxDZeta(indxlc) dxDZeta(indc) dxDZeta(indxcr) dxDbdy(5)];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
     
        HL=Dwl(1)+etaL;DZL=Dwl(1)+ZshipL;
        HR=Dwl(2)+etaR;DZR=Dwl(2)+ZshipR;
        
        [AmatD,BvectD] =fun_construct_matrices_diffraction_Ndom(Ndom,...
                0,0,tauGbdy,dxtauGbdy,Dbdy,dxDbdy,DZetabdy,dxDZetabdy,dr,...
                Evmodes,kappaS,etaWl,dxetaWl,phiWl,dxphiWl,Xbdy);
%             
%         [AmatD0,BvectD0] =fun_construct_matrices_diffraction_Ndom0(Ndom,...
%                 0,0,tauGbdy,dxtauGbdy,Dbdy,dxDbdy,DZetabdy,dxDZetabdy,dr,...
%                 Evmodes,kappaS,etaWl,dxetaWl,phiWl,dxphiWl,Xbdy);
            
       % [Amat, Bvect] =fun_construct_matrices_diffraction1(tauAzWl,taugWl,dxtauAzWl,dxtaugWl,Dwl,dxDwl,dr,Evmodes,kappaS,X0,phiWl,dxphiWl,etaWl,dxetaWl,Xwl);

        
        Coef=(AmatD\BvectD);
%         assignin('base','Coef',Coef);
%         assignin('base','Amat',AmatD);
%         assignin('base','Bvect',BvectD);
%         Coef0=AmatD0\BvectD0;
%         assignin('base','Coef0',Coef0);
%         assignin('base','Amat0',AmatD0);
%         assignin('base','Bvect0',BvectD0);
       
 
        
        Fev=zeros(size(phi));FevL=0;FevR=0;
        Lev=zeros(size(phi));
        
        tau=zeros(size(chiSWl));
        tauZ0=zeros(size(chiSWl));
        Nm=Ndom*2-2;
        idAd=0;idCoef1=2;idC=Ndom/2;idAdE=0;idCoef2=Nm-2;
       
     
        InQ0TEtaL=fun_IntegQn(kappaS(1,1),etaL,Dbdy(1),-dr(1),etaL);
        InQ0TEtaR=fun_IntegQn(kappaS(1,2),etaR,Dbdy(end),-dr(end),etaR);
       
        betaxd=phiL*InQ0TEtaL-phiR*InQ0TEtaR;
        
        InQ0TEtaL=fun_IntegQn(kappaS(1,1),etaL,Dbdy(1),-dr(1),etaL);
        InQ0TEtaR=fun_IntegQn(kappaS(1,2),etaR,Dbdy(end),-dr(end),etaR);
        InZQ0TEtaL=fun_IntegZQn(kappaS(1,1),etaL,Dbdy(1),-dr(1),etaL);
        InZQ0TEtaR=fun_IntegZQn(kappaS(1,2),etaR,Dbdy(end),-dr(end),etaR);
        
        betathetad=-phiL.*(InZQ0TEtaL-Zc*InQ0TEtaL)+phiR.*(InZQ0TEtaR-Zc*InQ0TEtaR);
        
        for jj=1:(Ndom-2)/2
            if jj==(Ndom-2)/2, idAdE=1; end
            int1=indxBdy(jj)+idAd:indxBdy(jj+1);
            tau(int1)=tau(int1)+Coef(idCoef1).*taug(int1)+Coef(idCoef1-1);
            tauZ0(int1)=tauZ0(int1)+Coef(idCoef1).*taug(int1)+Coef(idCoef1-1);
            int2=indxBdy(end-jj)+idAdE:indxBdy(end-jj+1)-idAd;
            tau(int2)=tau(int2)+Coef(idCoef2).*taug(int2)+Coef(idCoef2-1);
            tauZ0(int2)=tauZ0(int2)+Coef(idCoef2).*taug(int2)+Coef(idCoef2-1);
            
            if jj<(Ndom-2)/2
                Intvarphi0Zetad_dprev=fun_IntegVarphin(0,Dbdy(jj+1),-dr(jj+1),-dr(jj));
                betaxd=betaxd+(Coef(idCoef1).*taug(indxBdy(jj+1))+Coef(idCoef1-1))*Intvarphi0Zetad_dprev;
                
                Intvarphi0Zetad_dnext=fun_IntegVarphin(0,Dbdy(end-jj),-dr(end-jj),-dr(end-jj+1));
                betaxd=betaxd-(Coef(idCoef2).*taug(indxBdy(end-jj))+Coef(idCoef2-1))*Intvarphi0Zetad_dnext;
                
                IntZvarphi0Zetad_dprev=fun_IntegZVarphin(0,Dbdy(jj+1),-dr(jj+1),-dr(jj));
                betathetad=betathetad-(Coef(idCoef1).*taug(indxBdy(jj+1))+Coef(idCoef1-1))*(IntZvarphi0Zetad_dprev-Zc*Intvarphi0Zetad_dprev);
                
                IntZvarphi0Zetad_dnext=fun_IntegZVarphin(0,Dbdy(end-jj),-dr(end-jj),-dr(end-jj+1));
                betathetad=betathetad+(Coef(idCoef2).*taug(indxBdy(end-jj))+Coef(idCoef2-1))*(IntZvarphi0Zetad_dnext-Zc*Intvarphi0Zetad_dnext);
            end
            idAd=1;idCoef1=idCoef1+2; idCoef2=idCoef2-2;
        end
        
        
        
        %%%for calculating d\xi K(\phi,\xi)
        QLz0=cos(kappaS(1,1).*DZL)./cos(kappaS(1,1).*HL);
        QRz0=cos(kappaS(1,2).*DZR)./cos(kappaS(1,2).*HR);
        dxQLz0=kappaS(1,1).*(-dxDwl(1).*sin(kappaS(1,1).*DZL)+(dxDwl(1)+dxetaL).*tan(kappaS(1,1).*HL).*cos(kappaS(1,1).*DZL))./cos(kappaS(1,1).*HL);
        dxQRz0=kappaS(1,2).*(-dxDwl(2).*sin(kappaS(1,2).*DZR)+(dxDwl(2)+dxetaR).*tan(kappaS(1,2).*HR).*cos(kappaS(1,2).*DZR))./cos(kappaS(1,2).*HR);
        
        dzQLz0=-kappaS(1,1).*sin(kappaS(1,1).*DZL)./cos(kappaS(1,1).*HL);
        dzQRz0=-kappaS(1,2).*sin(kappaS(1,2).*DZR)./cos(kappaS(1,2).*HR);
        
        
        dPhi1dx_zL=dxphiL.*QLz0+phiL.*dxQLz0;
        dPhi3dx_zR=dxphiR.*QRz0+phiR.*dxQRz0;
        dPhi1dz_zL=phiL.*dzQLz0;
        dPhi3dz_zR=phiR.*dzQRz0;
        %%%%%%%%%%%%%%%%%%%%%
        
        if Evmodes>0
            for m=1:Evmodes
                eps1_m=Coef(Nm*m+1-2);
                epsNdom_m=Coef(Nm*(m+1)-2);
               
                phif(1:indWl(1))=phif(1:indWl(1))+eps1_m*exp(kappaS(m+1,1).*(x(1:indWl(1))-Xwl(1)));
                phif(indWl(2):end)=phif(indWl(2):end)+epsNdom_m*exp(-kappaS(m+1,2)*(x(indWl(2):end)-Xwl(2)));

                mu_m=m*pi./(DZeta);
                tauev=zeros(size(x));
                tauevZ0=zeros(size(x));
                %
                idAd=0;idAdE=0; idCoef1=Nm*m+2-2; idCoef2=Nm*m+Nm-2-2;
                InQmTEtaL=fun_IntegQn(kappaS(m+1,1),etaL,Dbdy(1),-dr(1),etaL);
                InQmTEtaR=fun_IntegQn(kappaS(m+1,2),etaR,Dbdy(end),-dr(end),etaR);
                betaxd=betaxd+eps1_m*InQmTEtaL-epsNdom_m*InQmTEtaR;
               
                InZQmTEtaL=fun_IntegZQn(kappaS(m+1,1),etaL,Dbdy(1),-dr(1),etaL);
                InZQmTEtaR=fun_IntegZQn(kappaS(m+1,2),etaR,Dbdy(end),-dr(end),etaR);
                betathetad=betathetad-eps1_m*(InZQmTEtaL-Zc*InQmTEtaL)+epsNdom_m*(InZQmTEtaR-Zc*InQmTEtaR);
                
                 
                 for jj=1:(Ndom-2)/2
                     if jj==(Ndom-2)/2, idAdE=1; end
                     int1=indxBdy(jj)+idAd:indxBdy(jj+1);
                     tauev(int1)=((Coef(idCoef1).*exp(mu_m(int1).*(x(int1)-Xbdy(jj)))...
                         +Coef(idCoef1+1).*exp(-mu_m(int1).*(x(int1)-Xbdy(jj+1)))).*cos(mu_m(int1).*(DZeta(int1))));
                     tauevZ0(int1)=((Coef(idCoef1).*exp(mu_m(int1).*(x(int1)-Xbdy(jj)))...
                         +Coef(idCoef1+1).*exp(-mu_m(int1).*(x(int1)-Xbdy(jj+1)))).*cos(mu_m(int1).*(D(int1))));
                     int2=indxBdy(end-jj)+idAdE:indxBdy(end-jj+1)-idAd;
                     tauev(int2)=((Coef(idCoef2).*exp(mu_m(int2).*(x(int2)-Xbdy(end-jj)))...
                         +Coef(idCoef2+1).*exp(-mu_m(int2).*(x(int2)-Xbdy(end-jj+1)))).*cos(mu_m(int2).*(DZeta(int2))));
                     tauevZ0(int2)=((Coef(idCoef2).*exp(mu_m(int2).*(x(int2)-Xbdy(end-jj)))...
                         +Coef(idCoef2+1).*exp(-mu_m(int2).*(x(int2)-Xbdy(end-jj+1)))).*cos(mu_m(int2).*(D(int2))));
                     
                     if jj<(Ndom-2)/2
                         mu_mLd=m*pi./(Dbdy(jj)-dr(jj));
                         mu_mRd=m*pi./(Dbdy(end-jj+1)-dr(end-jj+1));
                         IntvarphimZetad_dprev=fun_IntegVarphin(mu_mLd,Dbdy(jj+1),-dr(jj+1),-dr(jj));
                         betaxd=betaxd+(Coef(idCoef1).*exp(mu_mLd.*(Xbdy(jj+1)-Xbdy(jj)))...
                             +Coef(idCoef1+1)).*IntvarphimZetad_dprev;
                         IntvarphimZetad_dnext=fun_IntegVarphin(mu_mRd,Dbdy(end-jj),-dr(end-jj),-dr(end-jj+1));
                         betaxd=betaxd-(Coef(idCoef2)+Coef(idCoef2+1).*exp(-mu_mRd.*(Xbdy(end-jj)-Xbdy(end-jj+1))))*IntvarphimZetad_dnext;
                         
                         IntZvarphimZetad_dprev=fun_IntegZVarphin(mu_mLd,Dbdy(jj+1),-dr(jj+1),-dr(jj));
                         betathetad=betathetad-(Coef(idCoef1).*exp(mu_mLd.*(Xbdy(jj+1)-Xbdy(jj)))...
                             +Coef(idCoef1+1)).*(IntZvarphimZetad_dprev-Zc*IntvarphimZetad_dprev);
                         
                         IntZvarphimZetad_dnext=fun_IntegZVarphin(mu_mRd,Dbdy(end-jj),-dr(end-jj),-dr(end-jj+1));
                         
                         betathetad=betathetad+(Coef(idCoef2)+Coef(idCoef2+1).*...
                             exp(-mu_mRd.*(Xbdy(end-jj)-Xbdy(end-jj+1))))*(IntZvarphimZetad_dnext-Zc*IntvarphimZetad_dnext);
                         
                     end
                     
                     idAd=1;idCoef1=idCoef1+2;idCoef2=idCoef2-2;
                 end
                 tauev=tauev.*chiship(:,1);
                 tauev(indxBdy(1))=Coef(Nm*m+1-2);
                 tauev(indxBdy(end))=Coef(Nm*(m+1)-2);
                 tau=tau+tauev;
                
                 tauevZ0=tauevZ0.*chiship(:,1);
                 tauevZ0(indxBdy(1))=Coef(Nm*m+1-2);
                 tauevZ0(indxBdy(end))=Coef(Nm*(m+1)-2);
                 tauZ0=tauZ0+tauevZ0;
                 
                 
                 
                 
                IntQmDetaL=fun_IntegQn(kappaS(m+1,1),etaWl(1),Dwl(1),-Dwl(1),etaWl(1));
                IntQmDetaR=fun_IntegQn(kappaS(m+1,2),etaWl(2),Dwl(2),-Dwl(2),etaWl(2));
                IntdxQmDetaL=fun_IntegdxQn(kappaS(m+1,1),etaWl(1),Dwl(1),dxetaWl(1),dxDwl(1),-Dwl(1),etaWl(1));
                IntdxQmDetaR=fun_IntegdxQn(kappaS(m+1,2),etaWl(2),Dwl(2),dxetaWl(2),dxDwl(2),-Dwl(2),etaWl(2));

                
%                 dxIntQmDetaL_x=gradient(IntQmDetaL_x,dx);
%                 dxIntdxQmDetaL_x=gradient(IntdxQmDetaL_x,dx);
               
%                 dxIntQmDetaR_x=gradient(IntQmDetaR_x,dx);
%                 dxIntdxQmDetaR_x=gradient(IntdxQmDetaR_x,dx);
                
                FevL=FevL+eps1_m*(kappaS(m+1,1).*IntQmDetaL+IntdxQmDetaL);
                FevR=FevR+epsNdom_m*(-kappaS(m+1,2).*IntQmDetaR+IntdxQmDetaR);  
%                 
%                  Lev(1:indl)=Lev(1:indl)-eps1_m.*exp(kappaS(m+1,1).*(x(1:indl)-X0+b)).*(kappaS(m+1,1).^2.*IntQmDetaL_x(1:indl)+kappaS(m+1,1).*dxIntQmDetaL_x(1:indl)+kappaS(m+1,1).*IntdxQmDetaL_x(1:indl)+dxIntdxQmDetaL_x(1:indl));
%                  Lev(indr:end)=Lev(indr:end)-eps3_m.*exp(-kappaS(m+1,2).*(x(indr:end)-X0-b)).*(kappaS(m+1,2).^2.*IntQmDetaR_x(indr:end)-kappaS(m+1,2).*dxIntQmDetaR_x(indr:end)-kappaS(m+1,2).*IntdxQmDetaR_x(indr:end)+dxIntdxQmDetaR_x(indr:end));

%                 FevL=FevL+eps1_m*tan(kappaS(m+1,1)*HL)+eps1_m*(tan(kappaS(m+1,1)*HL).^2.*(dxDwl(1)+dxetaL)+kappaS(m+1,1).*dxDwl(1).*(1-1./cos(kappaS(m+1,1).*HL)));
%                 FevR=FevR-eps3_m*tan(kappaS(m+1,2)*HR)+eps3_m*(tan(kappaS(m+1,2)*HR).^2.*(dxDwl(2)+dxetaR)+kappaS(m+1,2).*dxDwl(2).*(1-1./cos(kappaS(m+1,2).*HR)));                 
%                 FevL=FevL+eps1_m*tan(kappaS(m+1,1)*HL);
%                 FevR=FevR-eps3_m*tan(kappaS(m+1,2)*HR);
                  Lev(1:indWl(1))=Lev(1:indWl(1))-eps1_m.*exp(kappaS(m+1,1).*(x(1:indWl(1))-Xwl(1))).*kappaS(m+1,1).*tan(kappaS(m+1,1).*D(1:indWl(1)));           
                  Lev(indWl(2):end)=Lev(indWl(2):end)-epsNdom_m.*exp(-kappaS(m+1,2).*(x(indWl(2):end)-Xwl(2))).*kappaS(m+1,2).*tan(kappaS(m+1,2).*D(indWl(2):end));   
               

                %%for calculating d\xi K(\phi,\xi)
                QLzm=cos(kappaS(m+1,1).*DZL)./cos(kappaS(m+1,1).*HL);
                QRzm=cos(kappaS(m+1,2).*DZR)./cos(kappaS(m+1,2).*HR);
                dxQLzm=kappaS(m+1,1).*(-dxDwl(1).*(sin(kappaS(m+1,1).*DZL))+(dxDwl(1)+dxetaL).*tan(kappaS(m+1,1).*HL).*cos(kappaS(m+1,1).*DZL))./cos(kappaS(m+1,1).*HL);
                dxQRzm=kappaS(m+1,2).*(-dxDwl(2).*(sin(kappaS(m+1,2).*DZR))+(dxDwl(2)+dxetaR).*tan(kappaS(m+1,2).*HR).*cos(kappaS(m+1,2).*DZR))./cos(kappaS(m+1,2).*HR);
                dzQLzm=-kappaS(m+1,1).*sin(kappaS(m+1,1).*DZL)./cos(kappaS(m+1,1).*HL);
                dzQRzm=-kappaS(m+1,2).*sin(kappaS(m+1,2).*DZR)./cos(kappaS(m+1,2).*HR);
                
                dPhi1dx_zL=dPhi1dx_zL+eps1_m.*(kappaS(m+1,1).*QLzm+dxQLzm);
                dPhi3dx_zR=dPhi3dx_zR+epsNdom_m.*(-kappaS(m+1,2).*QRzm+dxQRzm);
                dPhi1dz_zL=dPhi1dz_zL+eps1_m.*dzQLzm;
                dPhi3dz_zR=dPhi3dz_zR+epsNdom_m.*dzQRzm;
                %%%%%%%%%%%%%%%%%%%%%
            end
        end
        dPhi2dx_x=gradient(tau,dx);%sig1/b.*ones(size(x));
        
        Fev(indl)=FevL;Fev(indr)=FevR;
        sig1C=Coef((Nm-2)/2);
        Fsi=(dxtauAz+sig1C*dxtaug).*(D-max(dr));
        Fs=Fs+real(Fsi).*(chiship(:,1)+chiwaterline(:,1));
%         InQ0DEtaL=fun_IntegQn(kappaS(1,1),etaL,Dbdy(1),-Dbdy(1),-dr(1));
%         InQ0DEtaR=fun_IntegQn(kappaS(1,2),etaR,Dbdy(end),-Dbdy(2),-dr(end));
       
        
%         FphiL=dxphiL.*(D(indl)-dr(1));%InQ0DEtaL;
%         FphiR=dxphiR.*(D(indr)-dr(end));%InQ0DEtaR;
%         FphiWL=zeros(size(phi));
       %FphiWL(indl)=FphiL;FphiWL(indr)=FphiR;
%         FphiWL=(dxtauAz+sig1C*dxtaug).*(D-0);
%         FphiWL=FphiWL.*(chiship(:,1)+chiwaterline(:,1));
       
%         [dxphiL Fs(indl) FphiWL(indl) Fs(indl)-FphiWL(indl)]
%         [dxphiR  Fs(indr) FphiWL(indr) Fs(indr)-FphiWL(indr)]
%         disp('---------')
        tau=real(tau);
        dPhi2dx_x=real(dPhi2dx_x);
        dPhi1dz_zL=real(dPhi1dz_zL);
        dPhi3dz_zR=real(dPhi3dz_zR);
        
        betathetad=betathetad+...
                trapz(x(indxBdy(1):indxBdy(end)),(x(indxBdy(1):indxBdy(end))-Xc).*tau(indxBdy(1):indxBdy(end)));
        indL=indxBdy(1);indR=indxBdy(end);
       
%         aa=((phi(indL)-tau(indL)).*x(indR)-(phi(indR)-tau(indR)).*x(indL))/(x(indR)-x(indL));
%         bb=((phi(indL)-tau(indL))-(phi(indR)-tau(indR)))/(x(indL)-x(indR));
%         tauS=tau+aa+bb*x;
         
        tauS=tau;%.*cosh(kappaS(1,1).*D)./cosh(kappaS(1,1).*DZeta);
        
        phif=phif.*(1-chiship(:,1))+tauS.*chiship(:,1);
        phifZ0=phif.*(1-chiship(:,1))+real(tauZ0).*chiship(:,1);
        beta.z(1)       = trapz(x(indl:indr),tau(indl:indr));
        beta.x(1)       = real(betaxd);
        beta.theta(1)   = real(betathetad);%trapz(x(indl:indr),phif(indl:indr).*nutild.theta(1,indl:indr)');
        dxi_K.z(1)      = trapz(x(indl:indr),dPhi2dx_x(indl:indr).^2)/2;%dPhi2dz is zeros (impermeability condition)
       
        dxi_K.x(1)      = trapz(ZshipL,(dPhi1dz_zL.^2)./2) ...
            -trapz(ZshipR,(dPhi3dz_zR.^2)./2); %% n3 is zero at x=-b and x=b and dphidx=0  (impermeability condition)
        dxi_K.theta(1)  = -trapz(ZshipL,(ZshipL-Zc).*(dPhi1dz_zL.^2)./2) ...
            +trapz(ZshipR,(ZshipR-Zc).*(dPhi3dz_zR.^2)./2)...
             +trapz(x(indl:indr),(x(indl:indr)-Xc).*(dPhi2dx_x(indl:indr).^2))/2;
     
%         dxi_K.x(1)      = (trapz(ZshipL,(dPhi1dx_zL.^2+dPhi1dz_zL.^2)./2-dPhi1dz_zL.*(dPhi1dx_zL.*dr./dx)./sqrt(1+(dr./dx).^2)) ...
%             -trapz(ZshipR,(dPhi3dx_zR.^2+dPhi3dz_zR.^2)./2-dPhi3dz_zR.*(-dPhi3dx_zR.*dr./dx)./sqrt(1+(dr./dx).^2))); %% n3 is zero at x=-b and x=b and dphidx=0 
%         dxi_K.theta(1)  = -trapz(ZshipL,(ZshipL-Zc).*((dPhi1dx_zL.^2+dPhi1dz_zL.^2)./2-dPhi1dz_zL.*(dPhi1dx_zL.*dr./dx)./sqrt(1+(dr./dx).^2))) ...
%             +trapz(ZshipR,(ZshipR-Zc).*((dPhi3dx_zR.^2+dPhi3dz_zR.^2)./2-dPhi3dz_zR.*(-dPhi3dx_zR.*dr./dx)./sqrt(1+(dr./dx).^2)))...
%             +trapz(x(indl:indr)-Xc,(x(indl:indr)-Xc-Xc).*(dPhi2dx_x(indl:indr).^2))/2;

    % if time>60
    %    figure(111);
    %    plot(x,phi,'r',x,phif,'b');
    %    xlim([-5 5]);
    %    plot_properties;
    % end

end


