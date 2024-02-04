function [phif,tau,Fs,Fwl,Lev,beta,dxi_K,etaWl,dxetaWl,tauAz]=funS_tangential_potential3_2barges(calclin,k,x,xShip,Sdraft,Slength,sZ,phi,eta,chiship,chiwaterline,depth,Evmodes,XZrot0,lambda_p,kappaS,taug,L0,F0,Nonlin,ChiAdj)
Nship=2;

%%%% find anzatz tangential potential (before apply impermeability condition)
chiSWl=(chiship(:,end)+chiwaterline(:,end));
% tauAz=fun_tangential_potential_anzats(k,x,L0,F0,taug,eta,phi,Nonlin,ChiAdj,xShip,chiSWl,Nship);
tauAz=zeros(size(phi));
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

ChiSWl1=(chiship(:,1)+chiwaterline(:,1));
ChiSWl2=(chiship(:,2)+chiwaterline(:,2));
Xwl=[xShip(1,1) xShip(1,3) xShip(2,1) xShip(2,3)];
indWl=closest(x,Xwl);
Dwl=D(indWl);
dxDwl=dxD(indWl);

dr=Sdraft-sZ.';%draft

Z0=XZrot0(:,2);%%
X0=[xShip(1,2) xShip(2,2)];%%
kappaS(1,:)=-kappaS(1,:);

        %%%%%%%%%%%% if ship not moving in horizontal direction%%%%%%%%%%%
 
        if  lambda_p/sum(Slength)<10
            phi1L=(phi(indWl(1)-4)+phi(indWl(1)-3)+phi(indWl(1)-2)+phi(indWl(1)-1))/4;%
            phi1R=(phi(indWl(2)+1)+phi(indWl(2)+2))/2;%
            phi2L=(phi(indWl(3)-2)+phi(indWl(3)-1))/2;%
            phi2R=(phi(indWl(4)+1)+phi(indWl(4)+2)+phi(indWl(4)+3)+phi(indWl(4)+4))/4;%
            
            eta1L=(eta(indWl(1)-4)+eta(indWl(1)-3)+eta(indWl(1)-2)+eta(indWl(1)-1))/4;%
            eta1R=(eta(indWl(2)+1)+eta(indWl(2)+2))/2;%
            eta2L=(eta(indWl(3)-2)+eta(indWl(3)-1))/2;%
            eta2R=(eta(indWl(4)+1)+eta(indWl(4)+2)+eta(indWl(4)+3)+eta(indWl(4)+4))/4;%
            
            dxphi1L=(dxphi(indWl(1)-4)+dxphi(indWl(1)-3)+dxphi(indWl(1)-2)+dxphi(indWl(1)-1))/4;%
            dxphi1R=(dxphi(indWl(2)+1)+dxphi(indWl(2)+2))/2;%
            dxphi2L=(dxphi(indWl(3)-2)+dxphi(indWl(3)-1))/2;%
            dxphi2R=(dxphi(indWl(4)+1)+dxphi(indWl(4)+2)+dxphi(indWl(4)+3)+dxphi(indWl(4)+4))/4;%
            
            dxeta1L=(dxeta(indWl(1)-4)+dxeta(indWl(1)-3)+dxeta(indWl(1)-2)+dxeta(indWl(1)-1))/4;%
            dxeta1R=(dxeta(indWl(2)+1)+dxeta(indWl(2)+2))/2;%
            dxeta2L=(dxeta(indWl(3)-2)+dxeta(indWl(3)-1))/2;%
            dxeta2R=(dxeta(indWl(4)+1)+dxeta(indWl(4)+2)+dxeta(indWl(4)+3)+dxeta(indWl(4)+4))/4;%
            
            phiWl=[phi1L phi1R phi2L phi2R];
            etaWl=[eta1L eta1R eta2L eta2R];
            dxphiWl=[dxphi1L dxphi1R dxphi2L dxphi2R];
            dxetaWl=[dxeta1L dxeta1R dxeta2L dxeta2R];
            
        else
            phiWl=phi(indWl);
            etaWl=eta(indWl);
            dxphiWl=dxphi(indWl);%
            dxetaWl=dxeta(indWl);
        end
        
        if calclin==1
           etaWl=[0 0 0 0];
           dxetaWl=[0 0 0 0];      
        end
        tauAzWl=tauAz(indWl);      
        taugWl=taug(indWl); 
        dxtauAz=gradient(tauAz,dx);
        dxtaug=gradient(taug,dx);        dxtauAzWl=dxtauAz(indWl);
        dxtaugWl=dxtaug(indWl);
        
        [Amat, Bvect] =fun_construct_matrices_diffraction_2barges(tauAzWl,taugWl,dxtauAzWl,dxtaugWl,Dwl,dxDwl,dr,Evmodes,kappaS,X0,phiWl,dxphiWl,etaWl,dxetaWl,Xwl);
               
        Coef=Amat\Bvect;
          
%         Zship1L=linspace(-dr(1),eta1L,100);
%         Zship1R=linspace(-dr(1),eta1R,100);
%         Zship2L=linspace(-dr(2),eta2L,100);
%         Zship2R=linspace(-dr(2),eta2R,100);
%         
%         H1L=Dwl(1)+eta1L;DZ1L=Dwl(1)+Zship1L;
%         H1R=Dwl(2)+eta1R;DZ1R=Dwl(2)+Zship1R;
%         H2L=Dwl(3)+eta2L;DZ2L=Dwl(3)+Zship2L;
%         H2R=Dwl(4)+eta2R;DZ2R=Dwl(4)+Zship2R;
       
%         if max(isnan(Coef))==1
%          tau0=0;tau1=0;tau2=0;tau3=0;
%         else
        sig0=Coef(1);sig1=Coef(2);
        sig2=Coef(3);sig3=Coef(4);
%         end
        
        tau=(tauAz+(sig0+sig1.*taug)).*ChiSWl1+(tauAz+(sig2+sig3.*taug)).*ChiSWl2;
        Fwl=zeros(size(phi));Fev1L=0;Fev1R=0;Fev2L=0;Fev2R=0;
        Lev=zeros(size(phi));
        
        IntQ0Teta1L=fun_IntegQn(kappaS(1,1),etaWl(1),Dwl(1),-dr(1),etaWl(1));
        IntQ0Teta1R=fun_IntegQn(kappaS(1,2),etaWl(2),Dwl(2),-dr(1),etaWl(2));
        IntQ0Teta2L=fun_IntegQn(kappaS(1,3),etaWl(3),Dwl(3),-dr(2),etaWl(3));
        IntQ0Teta2R=fun_IntegQn(kappaS(1,4),etaWl(4),Dwl(4),-dr(2),etaWl(4));

        beta1xm=phiWl(1)*IntQ0Teta1L-phiWl(2)*IntQ0Teta1R;
        beta2xm=phiWl(3)*IntQ0Teta2L-phiWl(4)*IntQ0Teta2R;
        

        IntQ0Deta1L=fun_IntegQn(kappaS(1,1),etaWl(1),Dwl(1),-Dwl(1),etaWl(1));
        IntQ0Deta1R=fun_IntegQn(kappaS(1,2),etaWl(2),Dwl(2),-Dwl(2),etaWl(2));
        IntQ0Deta2L=fun_IntegQn(kappaS(1,3),etaWl(3),Dwl(3),-Dwl(3),etaWl(3));
        IntQ0Deta2R=fun_IntegQn(kappaS(1,4),etaWl(4),Dwl(4),-Dwl(4),etaWl(4));
        
        IntdxQ0Deta1L=fun_IntegdxQn(kappaS(1,1),etaWl(1),Dwl(1),dxetaWl(1),dxDwl(1),-Dwl(1),etaWl(1));
        IntdxQ0Deta1R=fun_IntegdxQn(kappaS(1,2),etaWl(2),Dwl(2),dxetaWl(2),dxDwl(2),-Dwl(2),etaWl(2));
        IntdxQ0Deta2L=fun_IntegdxQn(kappaS(1,3),etaWl(3),Dwl(3),dxetaWl(3),dxDwl(3),-Dwl(3),etaWl(3));
        IntdxQ0Deta2R=fun_IntegdxQn(kappaS(1,4),etaWl(4),Dwl(4),dxetaWl(4),dxDwl(4),-Dwl(4),etaWl(4));
        
        Fwl1=dxphiWl(1)*IntQ0Deta1L+phiWl(1)*IntdxQ0Deta1L;
        Fwl2=dxphiWl(2)*IntQ0Deta1R+phiWl(2)*IntdxQ0Deta1R;
        Fwl3=dxphiWl(3)*IntQ0Deta2L+phiWl(3)*IntdxQ0Deta2L;
        Fwl4=dxphiWl(4)*IntQ0Deta2R+phiWl(4)*IntdxQ0Deta2R;
        
 
        if Evmodes>0
            for jj=1:Evmodes
                eps1_m=Coef(4*Nship*(jj-1)+5);eps21_m=Coef(4*Nship*(jj-1)+6);eps22_m=Coef(4*Nship*(jj-1)+7);eps31_m=Coef(4*Nship*(jj-1)+8);
                eps32_m=Coef(4*Nship*(jj-1)+9);eps41_m=Coef(4*Nship*(jj-1)+10);eps42_m=Coef(4*Nship*(jj-1)+11);eps5_m=Coef(4*Nship*(jj-1)+12);

                muS1L_m=jj*pi/(Dwl(1)-dr(1));
                muS1R_m=jj*pi/(Dwl(2)-dr(1));
                muS2L_m=jj*pi/(Dwl(3)-dr(2));
                muS2R_m=jj*pi/(Dwl(4)-dr(2));
                varphi_m=cos(jj*pi);% at z=-draft;
                
               
                phif(1:indWl(1))=phif(1:indWl(1))+eps1_m*exp(kappaS(jj+1,1)*(x(1:indWl(1))-Xwl(1)));
                phif(indWl(2):indWl(3))=(phif(indWl(2):indWl(3))+...
                    eps31_m*exp(-kappaS(jj+1,2)*(x(indWl(2):indWl(3))-Xwl(2)))...
                    +eps32_m*exp(kappaS(jj+1,3)*(x(indWl(2):indWl(3))-Xwl(3))));
              
                phif(indWl(4):end)=phif(indWl(4):end)+eps5_m*exp(-kappaS(jj+1,4)*(x(indWl(4):end)-Xwl(4)));
           
                tauEvm1     =(eps21_m.*exp(muS1L_m.*(x(indWl(1):indWl(2))-Xwl(1)))+...
                             eps22_m.*exp(-muS1R_m.*(x(indWl(1):indWl(2))-Xwl(2)))).*varphi_m;
                tauEvm2     =(eps41_m.*exp(muS2L_m.*(x(indWl(3):indWl(4))-Xwl(3)))+...
                             eps42_m.*exp(-muS2R_m.*(x(indWl(3):indWl(4))-Xwl(4)))).*varphi_m;
                       
                tau(indWl(1):indWl(2))=tau(indWl(1):indWl(2))+tauEvm1;
                tau(indWl(3):indWl(4))=tau(indWl(3):indWl(4))+tauEvm2;
                
                IntQmTeta1L=fun_IntegQn(kappaS(jj+1,1),etaWl(1),Dwl(1),-dr(1),etaWl(1));
                IntQmTeta1R=fun_IntegQn(kappaS(jj+1,2),etaWl(2),Dwl(2),-dr(1),etaWl(2));
                IntQmTeta2L=fun_IntegQn(kappaS(jj+1,3),etaWl(3),Dwl(3),-dr(2),etaWl(3));
                IntQmTeta2R=fun_IntegQn(kappaS(jj+1,4),etaWl(4),Dwl(4),-dr(2),etaWl(4));
                
                beta1xm=beta1xm+(eps1_m*IntQmTeta1L-(eps31_m+eps32_m*exp(kappaS(jj+1,2)*(Xwl(2)-Xwl(3)))).*IntQmTeta1R);
                beta2xm=beta2xm+((eps31_m*exp(-kappaS(jj+1,3)*(Xwl(3)-Xwl(2)))+eps32_m).*IntQmTeta2L-eps5_m.*IntQmTeta2R);
              
                IntQmDeta1L=fun_IntegQn(kappaS(jj+1,1),etaWl(1),Dwl(1),-Dwl(1),etaWl(1));
                IntQmDeta1R=fun_IntegQn(kappaS(jj+1,2),etaWl(2),Dwl(2),-Dwl(2),etaWl(2));
                IntQmDeta2L=fun_IntegQn(kappaS(jj+1,3),etaWl(3),Dwl(3),-Dwl(3),etaWl(3));
                IntQmDeta2R=fun_IntegQn(kappaS(jj+1,4),etaWl(4),Dwl(4),-Dwl(4),etaWl(4));
                
                IntdxQmDeta1L=fun_IntegdxQn(kappaS(jj+1,1),etaWl(1),Dwl(1),dxetaWl(1),dxDwl(1),-Dwl(1),etaWl(1));
                IntdxQmDeta1R=fun_IntegdxQn(kappaS(jj+1,2),etaWl(2),Dwl(2),dxetaWl(2),dxDwl(2),-Dwl(2),etaWl(2));
                IntdxQmDeta2L=fun_IntegdxQn(kappaS(jj+1,3),etaWl(3),Dwl(3),dxetaWl(3),dxDwl(3),-Dwl(3),etaWl(3));
                IntdxQmDeta2R=fun_IntegdxQn(kappaS(jj+1,4),etaWl(4),Dwl(4),dxetaWl(4),dxDwl(4),-Dwl(4),etaWl(4));
                 
                  IntQmDeta1_x=fun_IntegQn(kappaS(jj+1,1),eta(1:indWl(1)),D(1:indWl(1)),-D(1:indWl(1)),eta(1:indWl(1)));
                  IntdxQmDeta1_x=fun_IntegdxQn(kappaS(jj+1,1),eta(1:indWl(1)),D(1:indWl(1)),dxeta(1:indWl(1)),dxD(1:indWl(1)),-D(1:indWl(1)),eta(1:indWl(1)));
                  IntdxxQmDeta1_x=0;
                  IntQmDeta3_x=fun_IntegQn(kappaS(jj+1,2),eta(indWl(2):indWl(3)),D(indWl(2):indWl(3)),-D(indWl(2):indWl(3)),eta(indWl(2):indWl(3)));
                  IntdxQmDeta3_x=fun_IntegdxQn(kappaS(jj+1,2),eta(indWl(2):indWl(3)),D(indWl(2):indWl(3)),dxeta(indWl(2):indWl(3)),dxD(indWl(2):indWl(3)),-D(indWl(2):indWl(3)),eta(indWl(2):indWl(3)));
                  IntdxxQmDeta3_x=0;
                  IntQmDeta5_x=fun_IntegQn(kappaS(jj+1,4),eta(indWl(4):end),D(indWl(4):end),-D(indWl(4):end),eta(indWl(4):end));
                  IntdxQmDeta5_x=fun_IntegdxQn(kappaS(jj+1,4),eta(indWl(4):end),D(indWl(4):end),dxeta(indWl(4):end),dxD(indWl(4):end),-D(indWl(4):end),eta(indWl(4):end));
                  IntdxxQmDeta5_x=0;
%               
                
                Fev1L=Fev1L+eps1_m*(kappaS(jj+1,1).*IntQmDeta1L+IntdxQmDeta1L);
                Fev1R=Fev1R+eps31_m*(-kappaS(jj+1,2).*IntQmDeta1R+IntdxQmDeta1R)...
                           +eps32_m*exp(kappaS(jj+1,2)*(Xwl(2)-Xwl(3))).*(kappaS(jj+1,2).*IntQmDeta1R+IntdxQmDeta1R);
                Fev2L=Fev2L+eps31_m*exp(-kappaS(jj+1,3)*(Xwl(3)-Xwl(2))).*(-kappaS(jj+1,3).*IntQmDeta2L+IntdxQmDeta2L)...
                           +eps32_m*(kappaS(jj+1,3).*IntQmDeta2L+IntdxQmDeta2L); 
                Fev2R=Fev2R+eps5_m*(-kappaS(jj+1,4).*IntQmDeta2R+IntdxQmDeta2R);       

%                 
                 Qm_x_z_D=1./cos(kappaS(jj+1,1).*(D+eta));
                 Qm_x_z_eta1=1; Qm_x_z_eta3=1; Qm_x_z_eta5=1;
                
                 dxQm_x_z_D=tan(kappaS(jj+1,1)*(D+eta)).*kappaS(jj+1,1).*(dxD+dxeta)./cos(kappaS(jj+1,1)*(D+eta));
                 dxQm_x_z_eta=tan(kappaS(jj+1,1)*(D+eta)).*kappaS(jj+1,1).*(dxeta);
                 dxD1=dxD(1:indWl(1));
                 dxeta1=dxeta(1:indWl(1));
                 Qm_x_z_D1=Qm_x_z_D(1:indWl(1));
                 dxQm_x_z_D1=dxQm_x_z_D(1:indWl(1));
                 dxQm_x_z_eta1=dxQm_x_z_eta(1:indWl(1));
                 
                 dxD3=dxD(indWl(2):indWl(3));
                 dxeta3=dxeta(indWl(2):indWl(3));
                 Qm_x_z_D3=Qm_x_z_D(indWl(2):indWl(3));
                 dxQm_x_z_D3=dxQm_x_z_D(indWl(2):indWl(3));
                 dxQm_x_z_eta3=dxQm_x_z_eta(indWl(2):indWl(3));
                 
                 dxD5=dxD(indWl(4):end);
                 dxeta5=dxeta(indWl(4):end);
                 Qm_x_z_D5=Qm_x_z_D(indWl(4):end);
                 dxQm_x_z_D5=dxQm_x_z_D(indWl(4):end);
                 dxQm_x_z_eta5=dxQm_x_z_eta(indWl(4):end);
                
%                  Lev(1:indWl(1))=Lev(1:indWl(1))-...
%                                 (eps1_m*exp(kappaS(jj+1,1)*(x(1:indWl(1))-Xwl(1)))...
%                                 .*((kappaS(jj+1,1)*Qm_x_z_D1+dxQm_x_z_D1).*dxD1+(kappaS(jj+1,1)*Qm_x_z_eta1+dxQm_x_z_eta1).*dxeta1...
%                                 +(kappaS(jj+1,1).^2.*IntQmDeta1_x+2.*kappaS(jj+1,1).*IntdxQmDeta1_x+IntdxxQmDeta1_x)));
%                  Lev(indWl(2):indWl(3))=Lev(indWl(2):indWl(3))-... 
%                                  (eps31_m.*exp(-kappaS(jj+1,2)*(x(indWl(2):indWl(3))-Xwl(2)))...
%                                  .*((-kappaS(jj+1,2)*Qm_x_z_D3+dxQm_x_z_D3).*dxD3+(kappaS(jj+1,2)*Qm_x_z_eta3+dxQm_x_z_eta3).*dxeta3...
%                                    +(kappaS(jj+1,2).^2.*IntQmDeta3_x+2.*kappaS(jj+1,2).*IntdxQmDeta3_x+IntdxxQmDeta3_x))...
%                                  +eps32_m.*exp(kappaS(jj+1,2)*(x(indWl(2):indWl(3))-Xwl(3)))...
%                                  .*((kappaS(jj+1,2)*Qm_x_z_D3+dxQm_x_z_D3).*dxD3+(kappaS(jj+1,2)*Qm_x_z_eta3+dxQm_x_z_eta3).*dxeta3...
%                                    +(kappaS(jj+1,2).^2.*IntQmDeta3_x+2.*kappaS(jj+1,2).*IntdxQmDeta3_x+IntdxxQmDeta3_x)));
%                  Lev(indWl(4):end)=Lev(indWl(4):end)-...
%                                 (eps5_m*exp(-kappaS(jj+1,4)*(x(indWl(4):end)-Xwl(4)))...
%                                 .*((kappaS(jj+1,4)*Qm_x_z_D5+dxQm_x_z_D5).*dxD5+(kappaS(jj+1,4).*Qm_x_z_eta5+dxQm_x_z_eta5).*dxeta5...
%                                 +(kappaS(jj+1,4).^2.*IntQmDeta5_x+2.*kappaS(jj+1,4).*IntdxQmDeta5_x+IntdxxQmDeta5_x)));            
       
                  %%% linear & Flat bottom          
                 Lev(1:indWl(1))=Lev(1:indWl(1))-eps1_m.*exp(kappaS(jj+1,1).*(x(1:indWl(1))-Xwl(1))).*kappaS(jj+1,1).*tan(kappaS(jj+1,1).*D(1:indWl(1)));           
                 Lev(indWl(2):indWl(3))=Lev(indWl(2):indWl(3))-...
                     (eps31_m.*exp(-kappaS(jj+1,2).*(x(indWl(2):indWl(3))-Xwl(2)))+...
                       eps32_m.*exp(kappaS(jj+1,2).*(x(indWl(2):indWl(3))-Xwl(3)))).*kappaS(jj+1,2).*tan(kappaS(jj+1,2).*D(indWl(2):indWl(3)));
                  Lev(indWl(4):end)=Lev(indWl(4):end)-eps5_m.*exp(-kappaS(jj+1,4).*(x(indWl(4):end)-Xwl(4))).*kappaS(jj+1,4).*tan(kappaS(jj+1,4).*D(indWl(4):end));   
                 %%%%%%%%%%%%%%%%%%%%%
            end
        end

        Fwl(indWl(1))=0*Fwl1+Fev1L; Fwl(indWl(2))=0*Fwl2+Fev1R;
        Fwl(indWl(3))=0*Fwl3+Fev2L; Fwl(indWl(4))=0*Fwl4+Fev2R;
        
        Fsi1=(dxtauAz+sig1*dxtaug).*(D-dr(1));
        Fsi2=(dxtauAz+sig3*dxtaug).*(D-dr(2));
        Fs=Fs+real(Fsi1.*ChiSWl1+Fsi2.*ChiSWl2);
        tau=real(tau);
     
        phif=phif.*(1-chiship(:,end))+tau.*chiship(:,end);
        beta1zm=trapz(x(indWl(1):indWl(2)),tau(indWl(1):indWl(2)));
        beta2zm=trapz(x(indWl(3):indWl(4)),tau(indWl(3):indWl(4)));
        
        beta.z       = real([beta1zm;beta2zm]);%trapz(x(indl:indr)-xc,tau(indl:indr));
        beta.x       = real([beta1xm;beta2xm]);
        beta.theta   = [0 0];
    end
    % if time>60
    %    figure(111);
    %    plot(x,phi,'r',x,phif,'b');
    %    xlim([-5 5]);
    %    plot_properties;
    % end


