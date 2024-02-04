function [radPot,addedMass,dampC,Lpsiev,Fpsi,FpsiEv,dxi_Krad]=fun_radiation_potential_anzats_Interp(g,k,x,Oprt,taug,zeta,dxzeta,IdMotion,Nonlin,nutild,xShip,XZrot0,chiSWl,Ns,Idplot,bathtype)

% if Nonlin<=2
%    Nonlin=2;
% end
dx=x(2)-x(1);
if strcmpi(IdMotion,'Free')
    Nmotion=3;
else
    Nmotion=1;
end
   
psiH=zeros(size(zeta));
psiS=zeros(size(zeta));
psiP=zeros(size(zeta));
FpsiS=zeros(size(x));
FpsiH=zeros(size(x));
FpsiP=zeros(size(x)); 
  
Xc=xShip(1,2); Zc=XZrot0(1,2);
       
for jj=1:Nmotion
    if strcmpi(IdMotion,'Heave') || (strcmpi(IdMotion,'Free')&&jj==1) 
        mu=nutild.z(end,:).';
    elseif strcmpi(IdMotion,'Surge') || (strcmpi(IdMotion,'Free')&&jj==2) 
        mu=nutild.x(end,:).';
    elseif strcmpi(IdMotion,'Pitch')  || (strcmpi(IdMotion,'Free')&&jj==3) 
       mu=nutild.theta(end,:).'; %%%???
    end
    mubar=mu;
   % size(mu)
%     figure;
%     plot(x,mu,'r',x,zeta,'b',x,chiSWl,'b')
ShiftDelta1=0;
ShiftDelta2=0;
 if jj==1
     ShiftDelta1=1;
 elseif jj==2
     ShiftDelta1=1;
     ShiftDelta2=1;
end


    
    for i=1:Ns
        
        SL=xShip(i,3)-xShip(i,1);
        indxL=closest(x,xShip(i,1))-ShiftDelta1;
        indxR=closest(x,xShip(i,3))+ShiftDelta1;
        if Ns==2
            if i==1
                indxL=closest(x,xShip(i,1));
                indxR=closest(x,xShip(i,3));
            else
                indxL=closest(x,xShip(i,1));
                indxR=closest(x,xShip(i,3));
            end
        end
        
        taugL    =   taug(indxL);
        taugR    =   taug(indxR);
        %apply solvability condition <mubar,\chiS>=<mubar,tau>=0
        fl         = (1/(taugL-taugR))*trapz(x(indxL:indxR),chiSWl(indxL:indxR).*mu(indxL:indxR).*(taug(indxL:indxR)-taugR));
        fr         = (1/(taugL-taugR))*trapz(x(indxL:indxR),chiSWl(indxL:indxR).*mu(indxL:indxR).*(taug(indxL:indxR)-taugL));
        DDR=zeros(length(x),1);DDL=zeros(length(x),1);
        DDR(indxR)  = 1/(dx);  DDL(indxL)  = 1/(dx);
        mubar(indxL:indxR)=mubar(indxL:indxR)-(fl.*DDL(indxL:indxR)-fr.*DDR(indxL:indxR));
    end
    
   
   
    psibar=Ifft(funOprt_Linv_interp(g,k,Oprt,mubar,bathtype));
    
    psi=zeros(size(psibar));
    for i=1:Ns
        indxL=closest(x,xShip(i,1))-ShiftDelta2;
        indxR=closest(x,xShip(i,3))+ShiftDelta2;
        A=(psibar(indxR).*taug(indxL)-psibar(indxL).*taug(indxR))/(taug(indxR)-taug(indxL));
        B=(psibar(indxR)-psibar(indxL))/(taug(indxL)-taug(indxR));
        psi(indxL:indxR)=psibar(indxL:indxR)+A+B.*taug(indxL:indxR);
    end
    psi=psi.*chiSWl;
     if strcmpi(IdMotion,'Heave')  || (strcmpi(IdMotion,'Free')&&jj==1) 
         FpsiH(indxL+ShiftDelta2)=fl;FpsiH(indxR-ShiftDelta2)=fr;
         Ma33=trapz(x,psi);
         if Ma33<0
            psi=-psi; 
         end
         psiH=psi;
     elseif strcmpi(IdMotion,'Surge') || (strcmpi(IdMotion,'Free')&&jj==2)
         FpsiS(indxL+ShiftDelta2)=fl;FpsiS(indxR-ShiftDelta2)=fr;
         Ma11=trapz(x,-psi.*dxzeta);
         if Ma11<0
            psi=-psi; 
         end
         psiS=psi;
     elseif strcmpi(IdMotion,'Pitch')  || (strcmpi(IdMotion,'Free')&&jj==3)
         FpsiP(indxL+ShiftDelta2)=fl;FpsiP(indxR-ShiftDelta2)=fr;
         Ma55=trapz(x,psi.*((x-Xc)+(zeta-Zc).*dxzeta));
         if Ma55<0
             psi=-psi;
         end
         psiP=psi;
     end
     if Idplot==1
         figure;
         subplot(3,1,1);
         plot(x,mubar,'b',x,mu,'r',x,zeta,'k')
         xlabel('x [m]');legend({'$\bar{\mu}$','$\mu$',},'interpreter','latex')
         subplot(3,1,2);
         plot(x,psibar,'b',x,psi,'r')
         xlabel('x [m]');legend({'$\bar{\psi}$','$\psi$'},'interpreter','latex')
         subplot(3,1,3)
         Lpsi=Ifft(funOprt_L_interp(g,k,Oprt,zeta,psi,Nonlin,bathtype));
         plot(x,mu,'b',x,Lpsi,'r')
         xlabel('x [m]');legend('\mu','L\psi')
     end
end
radPot.surge=psiS;
radPot.heave=psiH;
radPot.pitch=psiP;

Ma11=trapz(x,-psiS.*dxzeta);
Ma33=trapz(x,psiH);
Ma55=trapz(x,psiP.*((x-Xc)+(zeta-Zc).*dxzeta));
addedMass.z=Ma33;
addedMass.x=Ma11;
addedMass.theta=Ma55;
  if strcmpi(IdMotion,'Free')
  Ma31=trapz(x,-psiH.*dxzeta); Ma51=trapz(x,-psiP.*dxzeta);
  Ma13=trapz(x,psiS); Ma53=trapz(x,psiP);
  Ma15=trapz(x,psiS.*((x-Xc)+(zeta-Zc).*dxzeta));
  Ma35=trapz(x,psiH.*((x-Xc)+(zeta-Zc).*dxzeta));
  
  
  addedMass.S(1).Mat=[Ma11 Ma13 Ma15;Ma31 Ma33 Ma35;Ma51 Ma53 Ma55];
  else
  addedMass.S(1).Mat=diag([Ma11 Ma33 Ma55]);
  end
  dampC.S(1).Mat=[];
  Fpsi.heave=FpsiH.';Fpsi.surge=FpsiS.';Fpsi.pitch=FpsiP.';
  FpsiEv.surge=0;FpsiEv.heave=0;FpsiEv.pitch=0;
  Lpsiev.surge=0;Lpsiev.heave=0; Lpsiev.pitch=0;
  % dxi_Krad.z(1)=0;dxi_Krad.x(1)=0;dxi_Krad.theta(1)=0;
   dxi_Krad.z(1)      = -(xShip(1,3)-xShip(1,1))/2;%-trapz(x(indl:indr),1)/2;%dPhi2dz is zeros (impermeability condition)
   dxi_Krad.x(1)      = trapz(x,-dxzeta)/2;
   nu5=((x-Xc)+(zeta-Zc).*dxzeta);absN=sqrt(1+dxzeta.^2);
   dxi_Krad.theta(1) = trapz(x,chiSWl.*nu5.*(((x-Xc).^2+(zeta-Zc).^2)-(x-Xc).*nu5./absN))/2;
   
   
%   indxL=closest(x,xShip(1,1));
%   indxR=closest(x,xShip(1,3));
%   dPhiHdx_x=gradient(psiH,dx);%sig1/b.*ones(size(x));
%   
%    dzQLz0=-kappaS(1,1).*sin(kappaS(1,1).*DZL)./cos(kappaS(1,1).*HL);
%    dzQRz0=-kappaS(1,2).*sin(kappaS(1,2).*DZR)./cos(kappaS(1,2).*HR);
%       
%   dPhi1dz_zL=psiS(indxL).*dzQLz0;
%   dPhi3dz_zR=psiS(indxR).*dzQRz0;
%         
%   dxi_Krad.z(1)      = trapz(x(indl:indr),dPhiHdx_x(indl:indr).^2)/2;%dPhi2dz is zeros (impermeability condition)
%   dxi_Krad.x(1)      = trapz(ZshipL,(dPhi1dz_zL.^2)./2) ...
%       -trapz(ZshipR,(dPhi3dz_zR.^2)./2); %% n3 is zero at x=-b and x=b and dphidx=0  (impermeability condition)
%   dxi_Krad.theta(1)  = -trapz(ZshipL,(ZshipL-Zc).*(dPhi1dz_zL.^2)./2) ...
%       +trapz(ZshipR,(ZshipR-Zc).*(dPhi3dz_zR.^2)./2)...
%       +trapz(x(indl:indr),(x(indl:indr)-Xc).*(dPhi2dx_x(indl:indr).^2))/2;
