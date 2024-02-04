function [phif,tau,Fs,Fev,Lev,beta,dxi_K,etaWl,dxetaWl,tauAz]=funS_tangential_potential3_1(calclin,k,x,xShip,Sdraft,Slength,sZ,phi,eta,chiship,chiwaterline,depth,Evmodes,XZrot0,lambda_p,kappaS,taug,L0,F0,Nonlin,ChiAdj)
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

b=Slength/2;
dr=Sdraft-sZ.';%dr
Z0=XZrot0(1,2);%%
X0=xShip(1,2);%%
kappaS(1,:)=-kappaS(1,:);

indl=indWl(1);indr=indWl(2);
        
        if  lambda_p/Slength<10
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
        
        dxtauAzWl=dxtauAz(indWl);
        dxtaugWl=dxtaug(indWl);
        
        
        ZshipL=linspace(-dr,etaL,100);
        ZshipR=linspace(-dr,etaR,100);
        HL=Dwl(1)+etaL;DZL=Dwl(1)+ZshipL;
        HR=Dwl(2)+etaR;DZR=Dwl(2)+ZshipR;
        
        [Amat, Bvect] =fun_construct_matrices_diffraction1(tauAzWl,taugWl,dxtauAzWl,dxtaugWl,Dwl,dxDwl,dr,Evmodes,kappaS,X0,phiWl,dxphiWl,etaWl,dxetaWl,Xwl);

        
        Coef=Amat\Bvect;
        assignin('base','Coef',Coef);
        assignin('base','Amat',Amat);
        assignin('base','Bvect',Bvect);
       
        sig0=Coef(1);
        sig1=Coef(2);
 
        
        tau=(tauAz+(sig0+sig1.*taug)).*chiSWl;
        Fev=zeros(size(phi));FevL=0;FevR=0;
        Lev=zeros(size(phi));
        
        IntQ0TetaL=fun_IntegQn(kappaS(1,1),etaWl(1),Dwl(1),-dr(1),etaWl(1));
        IntQ0TetaR=fun_IntegQn(kappaS(1,2),etaWl(2),Dwl(2),-dr(1),etaWl(2));

        betaxm=phiL*IntQ0TetaL-phiR*IntQ0TetaR;
        
%         X0=pitchXZrot0(1);
%         Z0=pitchXZrot0(2);
%        
        IZQ0TetaL = fun_IntegZQn(kappaS(1,1),etaWl(1),Dwl(1),-dr(1),etaWl(1));
        IZQ0TetaR = fun_IntegZQn(kappaS(1,2),etaWl(2),Dwl(2),-dr(1),etaWl(2));

        betathetam1=phiR.*(IZQ0TetaR-Z0.*IntQ0TetaR)-phiL.*(IZQ0TetaL-Z0.*IntQ0TetaL);%
        betathetam2=trapz(x(indl:indr),tau(indl:indr).*(x(indl:indr)-X0));%(-sig0*(X0*2*b)+sig1*2*b^2/3);
        
        
        %%%for calculating d\xi K(\phi,\xi)
        QLz0=cos(kappaS(1,1).*DZL)./cos(kappaS(1,1).*HL);
        QRz0=cos(kappaS(1,2).*DZR)./cos(kappaS(1,2).*HR);
        dxQLz0=kappaS(1,1).*(-dxDwl(1).*sin(kappaS(1,1).*DZL)+(dxDwl(1)+dxetaL).*tan(kappaS(1,1).*HL).*cos(kappaS(1,1).*DZL))./cos(kappaS(1,1).*HL);
        dxQRz0=kappaS(1,2).*(-dxDwl(2).*sin(kappaS(1,2).*DZR)+(dxDwl(2)+dxetaR).*tan(kappaS(1,2).*HR).*cos(kappaS(1,2).*DZR))./cos(kappaS(1,2).*HR);
        
        dzQLz0=-kappaS(1,1).*sin(kappaS(1,1).*DZL)./cos(kappaS(1,1).*HL);
        dzQRz0=-kappaS(1,2).*sin(kappaS(1,2).*DZR)./cos(kappaS(1,2).*HR);
        
         dPhi2dx_x=gradient(tau,dx);%sig1/b.*ones(size(x));
        
        dPhi1dx_zL=dxphiL.*QLz0+phiL.*dxQLz0;
        dPhi3dx_zR=dxphiR.*QRz0+phiR.*dxQRz0;
        dPhi1dz_zL=phiL.*dzQLz0;
        dPhi3dz_zR=phiR.*dzQRz0;
        %%%%%%%%%%%%%%%%%%%%%
        if Evmodes>0
            for jj=1:Evmodes
                eps1_m=Coef(4*(jj-1)+3);
                eps21_m=Coef(4*(jj-1)+4);
                eps22_m=Coef(4*(jj-1)+5);
                eps3_m=Coef(4*(jj-1)+6);
%                 muS1L_m=jj*pi/(Dwl(1)-dr(1));
%                 muS1R_m=jj*pi/(Dwl(2)-dr(1));
%                 varphi_m=cos(jj*pi);
                mu_m=jj*pi./(D-dr);
                varphi_m=cos(mu_m.*(D-dr));
                 
           
                phif(1:indWl(1))=phif(1:indWl(1))+eps1_m*exp(kappaS(jj+1,1).*(x(1:indWl(1))-Xwl(1)));
                phif(indWl(2):end)=phif(indWl(2):end)+eps3_m*exp(-kappaS(jj+1,2)*(x(indWl(2):end)-Xwl(2)));

%                 tauEvm1     =(eps21_m.*exp(muS1L_m.*(x(indWl(1):indWl(2))-Xwl(1)))+...
%                              eps22_m.*exp(-muS1R_m.*(x(indWl(1):indWl(2))-Xwl(2)))).*varphi_m;
%                

                tauEvm1     =(eps21_m.*exp(mu_m(indWl(1):indWl(2)).*(x(indWl(1):indWl(2))-Xwl(1)))+...
                             eps22_m.*exp(-mu_m(indWl(1):indWl(2)).*(x(indWl(1):indWl(2))-Xwl(2)))).*varphi_m(indWl(1):indWl(2));
               
                tau(indWl(1):indWl(2))=tau(indWl(1):indWl(2))+tauEvm1;
                
              
                 
                IntQmDetaL=fun_IntegQn(kappaS(jj+1,1),etaWl(1),Dwl(1),-Dwl(1),etaWl(1));
                IntQmDetaR=fun_IntegQn(kappaS(jj+1,2),etaWl(2),Dwl(2),-Dwl(2),etaWl(2));
                IntdxQmDetaL=fun_IntegdxQn(kappaS(jj+1,1),etaWl(1),Dwl(1),dxetaWl(1),dxDwl(1),-Dwl(1),etaWl(1));
                IntdxQmDetaR=fun_IntegdxQn(kappaS(jj+1,2),etaWl(2),Dwl(2),dxetaWl(2),dxDwl(2),-Dwl(2),etaWl(2));

                
%                 dxIntQmDetaL_x=gradient(IntQmDetaL_x,dx);
%                 dxIntdxQmDetaL_x=gradient(IntdxQmDetaL_x,dx);
               
%                 dxIntQmDetaR_x=gradient(IntQmDetaR_x,dx);
%                 dxIntdxQmDetaR_x=gradient(IntdxQmDetaR_x,dx);
                
                FevL=FevL+eps1_m*(kappaS(jj+1,1).*IntQmDetaL+IntdxQmDetaL);
                FevR=FevR+eps3_m*(-kappaS(jj+1,2).*IntQmDetaR+IntdxQmDetaR);  

%                 
%                  Lev(1:indl)=Lev(1:indl)-eps1_m.*exp(kappaS(jj+1,1).*(x(1:indl)-X0+b)).*(kappaS(jj+1,1).^2.*IntQmDetaL_x(1:indl)+kappaS(jj+1,1).*dxIntQmDetaL_x(1:indl)+kappaS(jj+1,1).*IntdxQmDetaL_x(1:indl)+dxIntdxQmDetaL_x(1:indl));
%                  Lev(indr:end)=Lev(indr:end)-eps3_m.*exp(-kappaS(jj+1,2).*(x(indr:end)-X0-b)).*(kappaS(jj+1,2).^2.*IntQmDetaR_x(indr:end)-kappaS(jj+1,2).*dxIntQmDetaR_x(indr:end)-kappaS(jj+1,2).*IntdxQmDetaR_x(indr:end)+dxIntdxQmDetaR_x(indr:end));

%                 FevL=FevL+eps1_m*tan(kappaS(jj+1,1)*HL)+eps1_m*(tan(kappaS(jj+1,1)*HL).^2.*(dxDwl(1)+dxetaL)+kappaS(jj+1,1).*dxDwl(1).*(1-1./cos(kappaS(jj+1,1).*HL)));
%                 FevR=FevR-eps3_m*tan(kappaS(jj+1,2)*HR)+eps3_m*(tan(kappaS(jj+1,2)*HR).^2.*(dxDwl(2)+dxetaR)+kappaS(jj+1,2).*dxDwl(2).*(1-1./cos(kappaS(jj+1,2).*HR)));                 
%                 FevL=FevL+eps1_m*tan(kappaS(jj+1,1)*HL);
%                 FevR=FevR-eps3_m*tan(kappaS(jj+1,2)*HR);
                  Lev(1:indWl(1))=Lev(1:indWl(1))-eps1_m.*exp(kappaS(jj+1,1).*(x(1:indWl(1))-Xwl(1))).*kappaS(jj+1,1).*tan(kappaS(jj+1,1).*D(1:indWl(1)));           
                  Lev(indWl(2):end)=Lev(indWl(2):end)-eps3_m.*exp(-kappaS(jj+1,2).*(x(indWl(2):end)-Xwl(2))).*kappaS(jj+1,2).*tan(kappaS(jj+1,2).*D(indWl(2):end));   
                    
                IntQmTetaL=fun_IntegQn(kappaS(jj+1,1),etaWl(1),Dwl(1),-dr(1),etaWl(1));
                IntQmTetaR=fun_IntegQn(kappaS(jj+1,2),etaWl(2),Dwl(2),-dr(1),etaWl(2));
                betaxm=betaxm+(eps1_m.*IntQmTetaL-eps3_m.*IntQmTetaR);

                IZQmTetaL = fun_IntegZQn(kappaS(jj+1,1),etaWl(1),Dwl(1),-dr(1),etaWl(1));
                IZQmTetaR = fun_IntegZQn(kappaS(jj+1,2),etaWl(2),Dwl(2),-dr(1),etaWl(2));

                betathetam1=betathetam1+eps3_m.*(IZQmTetaR-Z0.*IntQmTetaR)-eps1_m.*(IZQmTetaL-Z0.*IntQmTetaL);%
                betathetam2=betathetam2+trapz(x(indl:indr),tauEvm1.*(x(indl:indr)-X0));

                %%for calculating d\xi K(\phi,\xi)
                QLzm=cos(kappaS(jj+1,1).*DZL)./cos(kappaS(jj+1,1).*HL);
                QRzm=cos(kappaS(jj+1,2).*DZR)./cos(kappaS(jj+1,2).*HR);
                dxQLzm=kappaS(jj+1,1).*(-dxDwl(1).*(sin(kappaS(jj+1,1).*DZL))+(dxDwl(1)+dxetaL).*tan(kappaS(jj+1,1).*HL).*cos(kappaS(jj+1,1).*DZL))./cos(kappaS(jj+1,1).*HL);
                dxQRzm=kappaS(jj+1,2).*(-dxDwl(2).*(sin(kappaS(jj+1,2).*DZR))+(dxDwl(2)+dxetaR).*tan(kappaS(jj+1,2).*HR).*cos(kappaS(jj+1,2).*DZR))./cos(kappaS(jj+1,2).*HR);
                dzQLzm=-kappaS(jj+1,1).*sin(kappaS(jj+1,1).*DZL)./cos(kappaS(jj+1,1).*HL);
                dzQRzm=-kappaS(jj+1,2).*sin(kappaS(jj+1,2).*DZR)./cos(kappaS(jj+1,2).*HR);
                
                dPhi2dx_x=dPhi2dx_x+...
                    mu_m.*(1-(x-X0).*dxD./(D-dr)).*(eps21_m.*exp(mu_m.*(x-X0))-eps22_m.*exp(-mu_m.*(x-X0))).*varphi_m;
                dPhi1dx_zL=dPhi1dx_zL+eps1_m.*(kappaS(jj+1,1).*QLzm+dxQLzm);
                dPhi3dx_zR=dPhi3dx_zR+eps3_m.*(-kappaS(jj+1,2).*QRzm+dxQRzm);
                dPhi1dz_zL=dPhi1dz_zL+eps1_m.*dzQLzm;
                dPhi3dz_zR=dPhi3dz_zR+eps3_m.*dzQRzm;
                %%%%%%%%%%%%%%%%%%%%%
            end
        end

        Fev(indl)=FevL;Fev(indr)=FevR;
        Fsi=(dxtauAz+sig1*dxtaug).*(D-dr);
        Fs=Fs+real(Fsi).*(chiship(:,1)+chiwaterline(:,1));
        tau=real(tau);
        dPhi2dx_x=real(dPhi2dx_x);
        dPhi1dz_zL=real(dPhi1dz_zL);
        dPhi3dz_zR=real(dPhi3dz_zR);
        
        phif=phif.*(1-chiship(:,1))+tau.*chiship(:,1);
        beta.z(1)       = trapz(x(indl:indr),tau(indl:indr));
        beta.x(1)       = real(betaxm);
        beta.theta(1)   = real(betathetam1+betathetam2);%trapz(x(indl:indr),phif(indl:indr).*nutild.theta(1,indl:indr)');
        dxi_K.z(1)      = trapz(x(indl:indr),dPhi2dx_x(indl:indr).^2)/2;%dPhi2dz is zeros (impermeability condition)
        dxi_K.x(1)      = trapz(ZshipL,(dPhi1dz_zL.^2)./2) ...
            -trapz(ZshipR,(dPhi3dz_zR.^2)./2); %% n3 is zero at x=-b and x=b and dphidx=0  (impermeability condition)
        dxi_K.theta(1)  = -trapz(ZshipL,(ZshipL-Z0).*(dPhi1dz_zL.^2)./2) ...
            +trapz(ZshipR,(ZshipR-Z0).*(dPhi3dz_zR.^2)./2)...
            +trapz(x(indl:indr),(x(indl:indr)-X0).*(dPhi2dx_x(indl:indr).^2))/2;
        
%         dxi_K.x(1)      = (trapz(ZshipL,(dPhi1dx_zL.^2+dPhi1dz_zL.^2)./2-dPhi1dz_zL.*(dPhi1dx_zL.*dr./dx)./sqrt(1+(dr./dx)^2)) ...
%             -trapz(ZshipR,(dPhi3dx_zR.^2+dPhi3dz_zR.^2)./2-dPhi3dz_zR.*(-dPhi3dx_zR.*dr./dx)./sqrt(1+(dr./dx)^2))); %% n3 is zero at x=-b and x=b and dphidx=0 
%         dxi_K.theta(1)  = -trapz(ZshipL,(ZshipL-Z0).*((dPhi1dx_zL.^2+dPhi1dz_zL.^2)./2-dPhi1dz_zL.*(dPhi1dx_zL.*dr./dx)./sqrt(1+(dr./dx)^2))) ...
%             +trapz(ZshipR,(ZshipR-Z0).*((dPhi3dx_zR.^2+dPhi3dz_zR.^2)./2-dPhi3dz_zR.*(-dPhi3dx_zR.*dr./dx)./sqrt(1+(dr./dx)^2)))...
%             +trapz(x(indl:indr)-X0,(x(indl:indr)-X0-X0).*(dPhi2dx_x(indl:indr).^2))/2;

    % if time>60
    %    figure(111);
    %    plot(x,phi,'r',x,phif,'b');
    %    xlim([-5 5]);
    %    plot_properties;
    % end

end


