function [phif,Fs,Fev,Lev,beta,dxi_K,etaL,etaR,dxetaL,dxetaR]=funS_tangential_potential3(calclin,k,x,xShip,Sdraft,Slength,sZ,phi,eta,chiship,chiwaterline,depth,Nship,modes,KappaS,Amat,Bvect,nutild,pitchXZrot0)

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
D=depth.';
x=x.';
dxD=gradient(D,dx);
iterx=1;
    for ii=1:Nship
        %%%%%%%%%%%% if ship not moving in horizontal direction%%%%%%%%%%%
        KappaL=KappaS(:,iterx);
        KappaR=KappaS(:,iterx+1);
        iterx=iterx+2;
   
        kappa0L=KappaL(1);kappaL=KappaL(2:end);
        kappa0R=KappaS(1);kappaR=KappaS(2:end);
        lambda_p=2*pi./(kappa0L./(-1i));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        b=Slength(ii)/2;
        
        indl=closest(x,xShip(ii,1));indr=closest(x,xShip(ii,3));
        indc=closest(x,xShip(ii,2));
        
        xc=x(indc);%xShip(ii,2);
        draft=Sdraft(ii)-sZ(ii);
        
        DL=depth(indl);
        DR=depth(indr);
        dxDL=dxD(indl);
        dxDR=dxD(indr);
        
        if  lambda_p/Slength(ii)<10
            phiL=(phi(indl-4)+phi(indl-3)+phi(indl-2)+phi(indl-1))/4;%
            phiR=(phi(indr+1)+phi(indr+2)+phi(indr+3)+phi(indr+4))/4;%
            etaL=(eta(indl-4)+eta(indl-3)+eta(indl-2)+eta(indl-1))/4;%
            etaR=(eta(indr+1)+eta(indr+2)+eta(indr+3)+eta(indr+4))/4;%
            dxphiL=(dxphi(indl-4)+dxphi(indl-3)+dxphi(indl-2)+dxphi(indl-1))/4;%
            dxphiR=(dxphi(indr+1)+dxphi(indr+2)+dxphi(indr+3)+dxphi(indr+4))/4;%
            dxetaL=(dxeta(indl-4)+dxeta(indl-3)+dxeta(indl-2)+dxeta(indl-1))/4;%
            dxetaR=(dxeta(indr+1)+dxeta(indr+2)+dxeta(indr+3)+dxeta(indr+4))/4;%
        else
            phiL=phi(indl);
            phiR=phi(indr);
            etaL=eta(indl);
            etaR=eta(indr);
            dxphiL=dxphi(indl);%
            dxphiR=dxphi(indr);%
            dxetaL=dxeta(indl);
            dxetaR=dxeta(indr);
        end
        if calclin==1
            etaL=0;etaR=0;dxetaL=0;dxetaR=0;
        end
        
        ZshipL=linspace(-draft,etaL,100);
        ZshipR=linspace(-draft,etaR,100);
        HL=DL+etaL;DZL=DL+ZshipL;
        HR=DR+etaR;DZR=DR+ZshipR;
        
        if calclin==0
            % if calculated nonlinearly, then define Amat and Bvect each iteration
            [Amat, Bvect] =fun_construct_matrices_diffraction(DL,DR,dxDL,dxDR,phiL,phiR,dxphiL,dxphiR,etaL,etaR,dxetaL,dxetaR,draft,b,modes,KappaL,KappaR);
        else
%             Bvect(1)=Bvect(1)*phiL;
%             Bvect(2)=Bvect(2)*phiR;
%             if modes>0
%             Bvect(3:4:end)=Bvect(3:4:end)*phiL;
%             Bvect(4:4:end)=Bvect(4:4:end)*phiR;
%             %                  Bvect(5:4:end)=Bvect(5:4:end)*dxphiL;
%             %                  Bvect(6:4:end)=Bvect(6:4:end)*dxphiR; 
%             end
            [Amat, Bvect] =fun_construct_matrices_diffraction(DL,DR,dxDL,dxDR,phiL,phiR,dxphiL,dxphiR,0,0,0,0,draft,b,modes,KappaL,KappaR);
        end
        Coef=Amat\Bvect;
        assignin('base','Coef',Coef);
        assignin('base','Amat',Amat);
        assignin('base','Bvect',Bvect);
       
        if max(isnan(Coef))==1
         tau0=0;
         tau1=0;
        else
        tau0=Coef(1);
        tau1=Coef(2);
        end
        
        tau=zeros(size(phi));
        tau(indl:indr)=(tau0+tau1*(x(indl:indr)-xc)/b);
        Fev=zeros(size(phi));FevL=0;FevR=0;
        Lev=zeros(size(phi));
        
        IntQ0TetaL=fun_IntegQn(kappa0L,etaL,DL,-draft,etaL);
        IntQ0TetaR=fun_IntegQn(kappa0R,etaR,DR,-draft,etaR);
        betaxm=phiL*IntQ0TetaL-phiR*IntQ0TetaR;
        
        X0=pitchXZrot0(1);
        Z0=pitchXZrot0(2);
       
        IZQ0TetaL = fun_IntegZQn(kappa0L,etaL,DL,-draft,etaL);
        IZQ0TetaR = fun_IntegZQn(kappa0R,etaR,DR,-draft,etaR);

        betathetam1=phiR.*(IZQ0TetaR-Z0.*IntQ0TetaR)-phiL.*(IZQ0TetaL-Z0.*IntQ0TetaL);%
        betathetam2=(-tau0*(X0*2*b)+tau1*2*b^2/3);
        
        
        %%%for calculating d\xi K(\phi,\xi)
        QLz0=cos(kappa0L.*DZL)./cos(kappa0L.*HL);
        QRz0=cos(kappa0R.*DZR)./cos(kappa0R.*HR);
        dxQLz0=kappa0L.*(-dxDL.*sin(kappa0L.*DZL)+(dxDL+dxetaL).*tan(kappa0L.*HL).*cos(kappa0L.*DZL))./cos(kappa0L.*HL);
        dxQRz0=kappa0R.*(-dxDR.*sin(kappa0R.*DZR)+(dxDR+dxetaR).*tan(kappa0R.*HR).*cos(kappa0R.*DZR))./cos(kappa0R.*HR);
        
        dzQLz0=-kappa0L.*sin(kappa0L.*DZL)./cos(kappa0L.*HL);
        dzQRz0=-kappa0R.*sin(kappa0R.*DZR)./cos(kappa0R.*HR);
        
        dPhi2dx_x=tau1/b.*ones(size(x));
        
        dPhi1dx_zL=dxphiL.*QLz0+phiL.*dxQLz0;
        dPhi3dx_zR=dxphiR.*QRz0+phiR.*dxQRz0;
        dPhi1dz_zL=phiL.*dzQLz0;
        dPhi3dz_zR=phiR.*dzQRz0;
        %%%%%%%%%%%%%%%%%%%%%
        if modes>0
            for jj=1:modes
                eps1_m=Coef(4*(jj-1)+3);
                eps21_m=Coef(4*(jj-1)+4);
                eps22_m=Coef(4*(jj-1)+5);
                eps3_m=Coef(4*(jj-1)+6);
                mu_m=jj*pi./(D-draft);
                varphi_m=cos(mu_m.*(D-draft));
                phif(1:indl)=phif(1:indl)+eps1_m*exp(kappaL(jj)*(x(1:indl)-xc+b));
                phif(indr:end)=(phif(indr:end)+eps3_m*exp(-kappaR(jj)*(x(indr:end)-xc-b)));
                tauEvm          =(eps21_m.*exp(mu_m(indl:indr).*(x(indl:indr)-xc))+eps22_m.*exp(-mu_m(indl:indr).*(x(indl:indr)-xc))).*varphi_m(indl:indr);
                tau(indl:indr)=tau(indl:indr)+tauEvm;
                
                IntQmDetaL=fun_IntegQn(kappaL(jj),etaL,DL,-DL,etaL);
                IntQmDetaR=fun_IntegQn(kappaR(jj),etaR,DR,-DR,etaR);
                IntdxQmDetaL=fun_IntegdxQn(kappaL(jj),etaL,DL,dxetaL,dxDL,-DL,etaL);
                IntdxQmDetaR=fun_IntegdxQn(kappaR(jj),etaR,DR,dxetaR,dxDR,-DR,etaR);
                 
                IntQmDetaL_x=fun_IntegQn(kappaL(jj),eta,D,-D,eta);
                dxIntQmDetaL_x=gradient(IntQmDetaL_x,dx);
                IntdxQmDetaL_x=fun_IntegdxQn(kappaL(jj),eta,D,dxeta,dxD,-D,eta);
                dxIntdxQmDetaL_x=gradient(IntdxQmDetaL_x,dx);
                IntQmDetaR_x=fun_IntegQn(kappaR(jj),eta,D,-D,eta);
                dxIntQmDetaR_x=gradient(IntQmDetaR_x,dx);
                IntdxQmDetaR_x=fun_IntegdxQn(kappaR(jj),eta,D,dxeta,dxD,-D,eta);
                dxIntdxQmDetaR_x=gradient(IntdxQmDetaR_x,dx);
                
                FevL=FevL+eps1_m*kappaL(jj).*IntQmDetaL+eps1_m*IntdxQmDetaL;
                FevR=FevR-eps3_m*kappaR(jj).*IntQmDetaR+eps3_m*IntdxQmDetaR;

%                 
%                  Lev(1:indl)=Lev(1:indl)-eps1_m.*exp(kappaL(jj).*(x(1:indl)-xc+b)).*(kappaL(jj).^2.*IntQmDetaL_x(1:indl)+kappaL(jj).*dxIntQmDetaL_x(1:indl)+kappaL(jj).*IntdxQmDetaL_x(1:indl)+dxIntdxQmDetaL_x(1:indl));
%                  Lev(indr:end)=Lev(indr:end)-eps3_m.*exp(-kappaR(jj).*(x(indr:end)-xc-b)).*(kappaR(jj).^2.*IntQmDetaR_x(indr:end)-kappaR(jj).*dxIntQmDetaR_x(indr:end)-kappaR(jj).*IntdxQmDetaR_x(indr:end)+dxIntdxQmDetaR_x(indr:end));

%                 FevL=FevL+eps1_m*tan(kappaL(jj)*HL)+eps1_m*(tan(kappaL(jj)*HL).^2.*(dxDL+dxetaL)+kappaL(jj).*dxDL.*(1-1./cos(kappaL(jj).*HL)));
%                 FevR=FevR-eps3_m*tan(kappaR(jj)*HR)+eps3_m*(tan(kappaR(jj)*HR).^2.*(dxDR+dxetaR)+kappaR(jj).*dxDR.*(1-1./cos(kappaR(jj).*HR)));                 
%                 FevL=FevL+eps1_m*tan(kappaL(jj)*HL);
%                 FevR=FevR-eps3_m*tan(kappaR(jj)*HR);

                 Lev(1:indl)=Lev(1:indl)-eps1_m.*exp(kappaL(jj).*(x(1:indl)-xc+b)).*kappaL(jj).*tan(kappaL(jj).*(D(1:indl)));
                 Lev(indr:end)=Lev(indr:end)-eps3_m*exp(-kappaR(jj).*(x(indr:end)-xc-b)).*kappaR(jj).*tan(kappaR(jj).*(D(indr:end)));
              
                IntQmTetaL=fun_IntegQn(kappaL(jj),etaL,DL,-draft,etaL);
                IntQmTetaR=fun_IntegQn(kappaR(jj),etaR,DR,-draft,etaR);
                betaxm=betaxm+(eps1_m.*IntQmTetaL-eps3_m.*IntQmTetaR);

                IZQmTetaL = fun_IntegZQn(kappaL(jj),etaL,DL,-draft,etaL);
                IZQmTetaR = fun_IntegZQn(kappaR(jj),etaR,DR,-draft,etaR);

                betathetam1=betathetam1+eps3_m.*(IZQmTetaR-Z0.*IntQmTetaR)-eps1_m.*(IZQmTetaL-Z0.*IntQmTetaL);%
                betathetam2=betathetam2+trapz(x(indl:indr)-xc,tauEvm.*(x(indl:indr)-xc-X0));

                %%%for calculating d\xi K(\phi,\xi)
                QLzm=cos(kappaL(jj).*DZL)./cos(kappaL(jj).*HL);
                QRzm=cos(kappaR(jj).*DZR)./cos(kappaR(jj).*HR);
                dxQLzm=kappaL(jj).*(-dxDL.*(sin(kappaL(jj).*DZL))+(dxDL+dxetaL).*tan(kappaL(jj).*HL).*cos(kappaL(jj).*DZL))./cos(kappaL(jj).*HL);
                dxQRzm=kappaR(jj).*(-dxDR.*(sin(kappaR(jj).*DZR))+(dxDR+dxetaR).*tan(kappaR(jj).*HR).*cos(kappaR(jj).*DZR))./cos(kappaR(jj).*HR);
                dzQLzm=-kappaL(jj).*sin(kappaL(jj).*DZL)./cos(kappaL(jj).*HL);
                dzQRzm=-kappaR(jj).*sin(kappaR(jj).*DZR)./cos(kappaR(jj).*HR);
                
                dPhi2dx_x=dPhi2dx_x+...
                    mu_m.*(1-(x-xc).*dxD./(D-draft)).*(eps21_m.*exp(mu_m.*(x-xc))-eps22_m.*exp(-mu_m.*(x-xc))).*varphi_m;
                dPhi1dx_zL=dPhi1dx_zL+eps1_m.*(kappaL(jj).*QLzm+dxQLzm);
                dPhi3dx_zR=dPhi3dx_zR+eps3_m.*(-kappaR(jj).*QRzm+dxQRzm);
                dPhi1dz_zL=dPhi1dz_zL+eps1_m.*dzQLzm;
                dPhi3dz_zR=dPhi3dz_zR+eps3_m.*dzQRzm;
                %%%%%%%%%%%%%%%%%%%%%
            end
        end

        Fev(indl)=FevL;Fev(indr)=FevR;
        Fsi=tau1.*(D-draft)/b;
        Fs=Fs+real(Fsi).*(chiship(:,ii)+chiwaterline(:,ii));
        tau=real(tau);
        dPhi2dx_x=real(dPhi2dx_x);
        dPhi1dz_zL=real(dPhi1dz_zL);
        dPhi3dz_zR=real(dPhi3dz_zR);
        
        phif=phif.*(1-chiship(:,ii))+tau.*chiship(:,ii);
        beta.z(ii)       = trapz(x(indl:indr)-xc,tau(indl:indr));
        beta.x(ii)       = real(betaxm);
        beta.theta(ii)   = real(betathetam1+betathetam2);%trapz(x(indl:indr),phif(indl:indr).*nutild.theta(ii,indl:indr)');
        dxi_K.z(ii)      = trapz(x(indl:indr)-xc,dPhi2dx_x(indl:indr).^2)/2;%dPhi2dz is zeros (impermeability condition)
        dxi_K.x(ii)      = trapz(ZshipL,(dPhi1dz_zL.^2)./2) ...
            -trapz(ZshipR,(dPhi3dz_zR.^2)./2); %% n3 is zero at x=-b and x=b and dphidx=0  (impermeability condition)
        dxi_K.theta(ii)  = -trapz(ZshipL,(ZshipL-Z0).*(dPhi1dz_zL.^2)./2) ...
            +trapz(ZshipR,(ZshipR-Z0).*(dPhi3dz_zR.^2)./2)...
            +trapz(x(indl:indr)-xc,(x(indl:indr)-X0-xc).*(dPhi2dx_x(indl:indr).^2))/2;
%         
%         dxi_K.x(ii)      = (trapz(ZshipL,(dPhi1dx_zL.^2+dPhi1dz_zL.^2)./2-dPhi1dz_zL.*(dPhi1dx_zL.*draft./dx)./sqrt(1+(draft./dx)^2)) ...
%             -trapz(ZshipR,(dPhi3dx_zR.^2+dPhi3dz_zR.^2)./2-dPhi3dz_zR.*(-dPhi3dx_zR.*draft./dx)./sqrt(1+(draft./dx)^2))); %% n3 is zero at x=-b and x=b and dphidx=0 
%         dxi_K.theta(ii)  = -trapz(ZshipL,(ZshipL-Z0).*((dPhi1dx_zL.^2+dPhi1dz_zL.^2)./2-dPhi1dz_zL.*(dPhi1dx_zL.*draft./dx)./sqrt(1+(draft./dx)^2))) ...
%             +trapz(ZshipR,(ZshipR-Z0).*((dPhi3dx_zR.^2+dPhi3dz_zR.^2)./2-dPhi3dz_zR.*(-dPhi3dx_zR.*draft./dx)./sqrt(1+(draft./dx)^2)))...
%             +trapz(x(indl:indr)-xc,(x(indl:indr)-X0-xc).*(dPhi2dx_x(indl:indr).^2))/2;
    end
    % if time>60
    %    figure(111);
    %    plot(x,phi,'r',x,phif,'b');
    %    xlim([-5 5]);
    %    plot_properties;
    % end

end


