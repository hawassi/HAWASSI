function [psi,Ma,dampC,Lpsiev,Fpsi,FpsiEv,FpsiS,dxi_Krad]=funS_radiation_potential_2barges(shippar,x,xShip,Sdraft,sZ,chiXcZ0,chiWL,nu,Smotion,etaWl,dxetaWl,D)

psiS=zeros(size(x));
psiH=zeros(size(x));
psiP=zeros(size(x));
LpsievS=zeros(size(psiS));
LpsievH=zeros(size(psiH));
LpsievP=zeros(size(psiP));
FpsiEvS=zeros(size(x));FevLS=0;FevRS=0;
FpsiS=zeros(size(x));
FpsiEvH=zeros(size(x));FevLH=0;FevRH=0;
FpsiH=zeros(size(x));
FpsiEvP=zeros(size(x));FevLP=0;FevRP=0;
FpsiP=zeros(size(x));

dx=x(2)-x(1);
    
Idmotion=Smotion(1);
Ns=2;% Nship

dxD=gradient(D,dx);
ChiSWl1=(shippar.form.chi(:,1)+shippar.form.chiWl(:,1)).';
ChiSWl2=(shippar.form.chi(:,2)+shippar.form.chiWl(:,2)).';
Dwl=shippar.Dwl;
dxDwl=shippar.dxDwl;
XWl=[xShip(1,1) xShip(1,3) xShip(2,1) xShip(2,3)];
% indl1=closest(x,XWl(1)); indr1=closest(x,XWl(2));
% indl2=closest(x,XWl(3)); indr2=closest(x,XWl(4));

dr=Sdraft-sZ.';%draft
Evmodes=shippar.Evmode;
Z0=shippar.form.XZrot0(:,2);%%
X0=[xShip(1,2) xShip(2,2)];%%%%
kappaS=shippar.kappa;
kappaS(1,:)=-kappaS(1,:);
[Amat, BvectS, BvectH, BvectP] =...
    fun_construct_matrices_radiation_2barges(Dwl,dxDwl,dr,shippar.Evmode,kappaS,Idmotion,X0,Z0,etaWl,dxetaWl,XWl);
assignin('base','Amat',Amat)
IntQ0Teta1L=fun_IntegQn(kappaS(1,1),etaWl(1),Dwl(1),-dr(1),etaWl(1));
IntQ0Teta1R=fun_IntegQn(kappaS(1,2),etaWl(2),Dwl(2),-dr(1),etaWl(2));
IntQ0Teta2L=fun_IntegQn(kappaS(1,3),etaWl(3),Dwl(3),-dr(2),etaWl(3));
IntQ0Teta2R=fun_IntegQn(kappaS(1,4),etaWl(4),Dwl(4),-dr(2),etaWl(4));


if strcmpi(Idmotion,'Surge')|| strcmpi(Idmotion,'Free')
    CoefS=Amat\BvectS;
    A1S=CoefS(1);A3LS=CoefS(2);A3RS=CoefS(3);A5S=CoefS(4);
    tau0S=CoefS(5);tau1S=CoefS(6);tau2S=CoefS(7);tau3S=CoefS(8);
    psiS=(tau0S+tau1S.*(x-X0(1))).*ChiSWl1;
    psiS=psiS+(tau2S+tau3S.*(x-X0(2))).*ChiSWl2;
    
    f11S1=A1S.*IntQ0Teta1L-A3LS*IntQ0Teta1R;
    f11S2=A3RS.*IntQ0Teta2L-A5S*IntQ0Teta2R;
    FpsiS1=tau1S.*(D-dr(1));
    FpsiS2=tau3S.*(D-dr(2));
    FpsiS=FpsiS1.*ChiSWl1+FpsiS2.*ChiSWl2;
end

if strcmpi(Idmotion,'Heave')|| strcmpi(Idmotion,'Free')
    CoefH=Amat\BvectH;
    A1H=CoefH(1);A3LH=CoefH(2);A3RH=CoefH(3);A5H=CoefH(4);
    tau0H=CoefH(5);tau1H=CoefH(6);tau2H=CoefH(7);tau3H=CoefH(8);
    psiH=(((D-dr(1)).^2-(x-X0(1)).^2)./2./(D-dr(1))+tau0H+tau1H.*(x-X0(1))).*ChiSWl1;
    psiH=psiH+(((D-dr(2)).^2-(x-X0(2)).^2)./2./(D-dr(2))+tau2H+tau3H.*(x-X0(2))).*ChiSWl2;
    
    f33S1=trapz(x,psiH.*ChiSWl1);
    f33S2=trapz(x,psiH.*ChiSWl2);
    FpsiH1=-(x-X0(1))+tau1H*(D-dr(1))+(dxD./(D-dr(1))).*((D.^2+dr(1)^2)/2-dr(1).*D+(-(D-dr(1)).^2/3+(x-X0(1)).^2)./2);
    FpsiH2=-(x-X0(2))+tau3H*(D-dr(2))+(dxD./(D-dr(2))).*((D.^2+dr(2)^2)/2-dr(2).*D+(-(D-dr(2)).^2/3+(x-X0(2)).^2)./2); 
    FpsiH=FpsiH1.*ChiSWl1+FpsiH2.*ChiSWl2;
end

if strcmpi(Idmotion,'Pitch')|| strcmpi(Idmotion,'Free')
    CoefP=Amat\BvectP;
    A1P=CoefP(1);A3LP=CoefP(2);A3RP=CoefP(3);A5P=CoefP(4);
    tau0P=CoefP(5);tau1P=CoefP(6);tau2P=CoefP(7);tau3P=CoefP(8);
    psiP=(((D-dr(1)).^2.*(x-X0(1))-(x-X0(1)).^3/3)./2./(D-dr(1))+tau0P+tau1P.*(x-X0(1))).*ChiSWl1;
    psiP=psiP+(((D-dr(2)).^2.*(x-X0(2))-(x-X0(2)).^3/3)./2./(D-dr(2))+tau2P+tau3P.*(x-X0(2))).*ChiSWl2;
    
    
    IZQ0Teta1L = fun_IntegZQn(kappaS(1,1),etaWl(1),Dwl(1),-dr(1),etaWl(1));
    IZQ0Teta1R = fun_IntegZQn(kappaS(1,2),etaWl(2),Dwl(2),-dr(1),etaWl(2));
    IZQ0Teta2L = fun_IntegZQn(kappaS(1,3),etaWl(3),Dwl(3),-dr(2),etaWl(3));
    IZQ0Teta2R = fun_IntegZQn(kappaS(1,4),etaWl(4),Dwl(4),-dr(2),etaWl(4));
    
    
    f551S1=A3LP.*(IZQ0Teta1R-Z0(1).*IntQ0Teta1R)-A1P.*(IZQ0Teta1L-Z0(1).*IntQ0Teta1L);%
    f552S1=trapz(x,(x-X0(1)).*psiP.*ChiSWl1);
    f551S2=A5P.*(IZQ0Teta2R-Z0(2).*IntQ0Teta2R)-A3RP.*(IZQ0Teta2L-Z0(2).*IntQ0Teta2L);%
    f552S2=trapz(x,(x-X0(1)).*psiP.*ChiSWl2);
    
    FpsiP1=(D-dr(1)).^2/6-(x-X0(1)).^2/2+tau1P*(D-dr(1))+(dxD./(D-dr(1))).*((x-X0(1)).*((D.^2+dr(1)^2)/2-dr(1).*D)+(-(D-dr(1)).^2/3.*(x-X0(1))+(x-X0(1)).^3/3)./2);
    FpsiP2=(D-dr(1)).^2/6-(x-X0(1)).^2/2+tau1P*(D-dr(1))+(dxD./(D-dr(1))).*((x-X0(1)).*((D.^2+dr(1)^2)/2-dr(1).*D)+(-(D-dr(1)).^2/3.*(x-X0(1))+(x-X0(1)).^3/3)./2);
    FpsiP=FpsiP1.*ChiSWl1+FpsiP2.*ChiSWl2;
end

if Evmodes>0
    for m=1:Evmodes
        muS1L_m=m*pi/(Dwl(1)-dr(1));
        muS1R_m=m*pi/(Dwl(2)-dr(1));
        muS2L_m=m*pi/(Dwl(3)-dr(2));
        muS2R_m=m*pi/(Dwl(4)-dr(2));
        
        IntQmTeta1L=fun_IntegQn(kappaS(m+1,1),etaWl(1),Dwl(1),-dr(1),etaWl(1));
        IntQmTeta1R=fun_IntegQn(kappaS(m+1,2),etaWl(2),Dwl(2),-dr(1),etaWl(2));
        IntQmTeta2L=fun_IntegQn(kappaS(m+1,3),etaWl(3),Dwl(3),-dr(2),etaWl(3));
        IntQmTeta2R=fun_IntegQn(kappaS(m+1,4),etaWl(4),Dwl(4),-dr(2),etaWl(4));
        
        if strcmpi(Idmotion,'Surge')|| strcmpi(Idmotion,'Free')
            epsS1m=CoefS(4*m*Ns+1);epsS21m=CoefS(4*m*Ns+2);epsS22m=CoefS(4*m*Ns+3);
            epsS31m=CoefS(4*m*Ns+4);epsS32m=CoefS(4*m*Ns+5);
            epsS41m=CoefS(4*m*Ns+6);epsS42m=CoefS(4*m*Ns+7);epsS5m=CoefS(4*m*Ns+8);
            psiSev=(epsS21m.*exp(muS1L_m.*((x-XWl(1))))+epsS22m.*exp(-muS1R_m.*(x-XWl(2)))).*cos(muS1L_m.*(Dwl(1)-dr(1))).*ChiSWl1;
            psiSev= psiSev+ (epsS41m.*exp(muS2L_m.*((x-XWl(3))))+epsS42m.*exp(-muS2R_m.*(x-XWl(4)))).*cos(muS2L_m.*(Dwl(3)-dr(2))).*ChiSWl2;
            psiS=psiS+psiSev;
            
            f11S1=f11S1+(epsS1m.*IntQmTeta1L-epsS31m.*IntQmTeta1R);
            f11S2=f11S2+(epsS32m.*IntQmTeta2L-epsS5m.*IntQmTeta2R);
            
            
        end
        
        if strcmpi(Idmotion,'Heave')|| strcmpi(Idmotion,'Free')
            epsH1m=CoefH(4*m*Ns+1);epsH21m=CoefH(4*m*Ns+2);epsH22m=CoefH(4*m*Ns+3);
            epsH31m=CoefH(4*m*Ns+4);epsH32m=CoefH(4*m*Ns+5);
            epsH41m=CoefH(4*m*Ns+6);epsH42m=CoefH(4*m*Ns+7);epsH5m=CoefH(4*m*Ns+8);
            psiHev=(epsH21m.*exp(muS1L_m.*((x-XWl(1))))+epsH22m.*exp(-muS1R_m.*(x-XWl(2)))).*cos(muS1L_m.*(Dwl(1)-dr(1))).*ChiSWl1;
            psiHev= psiHev+ (epsH41m.*exp(muS2L_m.*((x-XWl(3))))+epsH42m.*exp(-muS2R_m.*(x-XWl(4)))).*cos(muS2L_m.*(Dwl(3)-dr(2))).*ChiSWl2;
            psiH=psiH+psiHev;
            
            f33S1=f33S1+trapz(x,psiHev.*ChiSWl1);
            f33S2=f33S2+trapz(x,psiHev.*ChiSWl2);
            
            
        end
        
        if strcmpi(Idmotion,'Pitch')|| strcmpi(Idmotion,'Free')
            epsP1m=CoefP(4*m*Ns+1);epsP21m=CoefP(4*m*Ns+2);epsP22m=CoefP(4*m*Ns+3);
            epsP31m=CoefP(4*m*Ns+4);epsP32m=CoefP(4*m*Ns+5);
            epsP41m=CoefP(4*m*Ns+6);epsP42m=CoefP(4*m*Ns+7);epsP5m=CoefP(4*m*Ns+8);
            psiPev=(epsP21m.*exp(muS1L_m.*((x-XWl(1))))+epsP22m.*exp(-muS1R_m.*(x-XWl(2)))).*cos(muS1L_m.*(Dwl(1)-dr(1))).*ChiSWl1;
            psiPev= psiPev+ (epsP41m.*exp(muS2L_m.*((x-XWl(3))))+epsP42m.*exp(-muS2R_m.*(x-XWl(4)))).*cos(muS2L_m.*(Dwl(3)-dr(2))).*ChiSWl2;
            psiP=psiP+psiPev;
            
            IZQmTeta1L = fun_IntegZQn(kappaS(m+1,1),etaWl(1),Dwl(1),-dr(1),etaWl(1));
            IZQmTeta1R = fun_IntegZQn(kappaS(m+1,2),etaWl(2),Dwl(2),-dr(1),etaWl(2));
            IZQmTeta2L = fun_IntegZQn(kappaS(m+1,3),etaWl(3),Dwl(3),-dr(2),etaWl(3));
            IZQmTeta2R = fun_IntegZQn(kappaS(m+1,4),etaWl(4),Dwl(4),-dr(2),etaWl(4));
            
            f551S1=f551S1+epsP31m.*(IZQmTeta1R-Z0(1).*IntQmTeta1R)-epsP1m.*(IZQmTeta1L-Z0(1).*IntQmTeta1L);%
            f552S1=f552S1+trapz(x,(x-X0(1)).*psiPev.*ChiSWl1);
            f551S2=f551S2+epsP5m.*(IZQmTeta2R-Z0(2).*IntQmTeta2R)-epsP32m.*(IZQmTeta2L-Z0(2).*IntQmTeta2L);%
            f552S2=f552S2+trapz(x,(x-X0(2)).*psiPev.*ChiSWl2);
            
        end
    end
end


%figure;
if strcmpi(Idmotion,'Surge')
    Ma33S1=0;Ma55S1=0;
    Ma11S1=real(f11S1);
    B11S1=imag(f11S1);B33S1=0;B55S1=0;
    Ma33S2=0;Ma55S2=0;
    Ma11S2=real(f11S2);
    B11S2=imag(f11S2);B33S2=0;B55S2=0;
    psiS=real(psiS);
   % plot(x,psiS)
    
elseif strcmpi(Idmotion,'Heave')
    Ma33S1=real(f33S1);B33S1=imag(f33S1);
    Ma11S1=0;Ma55S1=0;B55S1=0;B11S1=0;
    Ma33S2=real(f33S2);B33S2=imag(f33S2);
    Ma11S2=0;Ma55S2=0;B55S2=0;B11S2=0;
    psiH=real(psiH);
   % plot(x,psiH)
    
elseif strcmpi(Idmotion,'Pitch')
    f55S1=f551S1+f552S1;
    Ma55S1=real(f55S1);B55S1=imag(f55S1);B33S1=0;B11S1=0;
    Ma11S1=0;Ma33S1=0;
    
    f55S2=f551S2+f552S2;
    Ma55S2=real(f55S2);B55S2=imag(f55S2);B33S2=0;B11S2=0;
    Ma11S2=0;Ma33S2=0;
    psiP=real(psiP);
  %  plot(x,psiP)
else
 %   plot(x,real(psiS),'r',x,real(psiH),'g',x,real(psiP),'b')
end
%hold on;

Ma.z(1)=Ma33S1;
Ma.x(1)=Ma11S1;
Ma.theta(1)=Ma55S1;
Ma.z(2)=Ma33S2;
Ma.x(2)=Ma11S2;
Ma.theta(2)=Ma55S2;

if strcmpi(Idmotion,'Free')
Ma13S1=0;Ma31S1=0;Ma51S1=0;
Ma15S1=0;Ma35S1=0;Ma53S1=0;
Ma13S2=0;Ma31S2=0;Ma51S2=0;% must be calculated
Ma15S2=0;Ma35S2=0;Ma53S2=0;

B13S1=0;B31S1=0;B51S1=0;
B15S1=0;B35S1=0;B53S1=0;
B13S2=0;B31S2=0;B51S2=0;% must be calculated
B15S2=0;B35S2=0;B53S2=0;
     
Ma.S(1).Mat=[Ma11S1 Ma13S1 Ma15S1;Ma31S1 Ma33S1 Ma35S1; Ma51S1 Ma53S1 Ma55S1];
Ma.S(2).Mat=[Ma11S2 Ma13S2 Ma15S2;Ma31S2 Ma33S2 Ma35S2; Ma51S2 Ma53S2 Ma55S2];
dampC.S(1).Mat=([B11S1 B13S1 B15S1; B31S1 B33S1 B35S1;B51S1 B53S1 B55S1]);
dampC.S(2).Mat=([B11S2 B13S2 B15S2; B31S2 B33S2 B35S2;B51S2 B53S2 B55S2]);
else
Ma.S(1).Mat=diag([Ma11S1 Ma33S1 Ma55S1]);
Ma.S(2).Mat=diag([Ma11S2 Ma33S2 Ma55S2]);
dampC.S(1).Mat=diag([B11S1 B33S1 B55S1]);
dampC.S(2).Mat=diag([B11S2 B33S2 B55S2]);
end

psi.surge=real(psiS).';psi.heave=real(psiH).';
psi.pitch=real(psiP).';
Fpsi.surge=FpsiS;Fpsi.heave=FpsiH;
Fpsi.pitch=FpsiP;
Lpsiev.surge=LpsievS;Lpsiev.heave=LpsievH;
Lpsiev.pitch=LpsievP;
FpsiEv.surge=FpsiEvS;FpsiEv.heave=FpsiEvH;
FpsiEv.pitch=FpsiEvP;

dxi_Krad.z      = [0;0];
dxi_Krad.x      = [0;0];
dxi_Krad.theta  = [0;0];

