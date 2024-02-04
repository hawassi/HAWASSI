function [Flinwave,Fnonlinwave,FlinwaveHalfLeft,FlinwaveHalfRight,FnonlinwaveHalfLeft,FnonlinwaveHalfRight]...
          =fun_forceS(g,k,x,eta,phif,dphidtS,Oprt,modelNonlin,bathtype,Nship,chiship,paramShip)

%%%For calculating force
phif_hat=fft(phif);
if strcmp(bathtype,'F')
    L0=Oprt.L0;
    dxphiS=Ifft(1i.*k.*phif_hat);
    L0phiS=Ifft(L0.*phif_hat);
    if modelNonlin==1
        deletaKS=0;
    else
        deletaKS=0.5*(dxphiS.^2-L0phiS.^2); %second order
        if  modelNonlin==3
            L0etaL0phiS=Ifft(L0.*fft(eta.*L0phiS));
            deletaKS=deletaKS+ L0phiS.*(eta.*Ifft(1i.*k.*fft(dxphiS))+L0etaL0phiS);
        end
    end
    
else
    dxphiS=Ifft(1i.*k.*phif_hat);
    L0phiS=Ifft(funOprtBathy_L0(g,k,Oprt,phif));
    if modelNonlin==1
        deletaKS=0;
    else
        deletaKS=0.5*(dxphiS.^2-L0phiS.^2); %second order
        if  modelNonlin==3
            L0etaL0phiS=Ifft(funOprtBathy_L0(g,k,Oprt,eta.*L0phiS));
            deletaKS=deletaKS+ L0phiS.*(eta.*Ifft(1i.*k.*fft(dxphiS))+L0etaL0phiS);
        end
    end
end
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Flinwave    =zeros(2*Nship,1);
Fnonlinwave =zeros(2*Nship,1);
FlinwaveHalfLeft=zeros(2*Nship,1);
FlinwaveHalfRight=zeros(2*Nship,1);
FnonlinwaveHalfLeft=zeros(2*Nship,1);
FnonlinwaveHalfRight=zeros(2*Nship,1);
for i=1:Nship
    nuz         = chiship(:,i);%+chiWlship(:,end);
    nux         = -paramShip.form.dxshape(:,i);
    Flinwave(i)   =trapz(x,dphidtS.*nuz);
    Fnonlinwave(i)=trapz(x,deletaKS.*nuz);
    Flinwave(Nship+i)   =trapz(x,dphidtS.*nux);
    Fnonlinwave(Nship+i)=trapz(x,deletaKS.*nux);
    IndxL=closest(x,paramShip.form.xShip(i,1));
    IndxC=closest(x,paramShip.form.xShip(i,2));
    IndxR=closest(x,paramShip.form.xShip(i,3));
    nu3LC        = nuz;nu3CR        = nuz;
    nu3LC(IndxC+1:IndxR)=0;nu3CR(IndxL:IndxC)=0;
    FlinwaveHalfLeft(i)=trapz(x,dphidtS.*nu3LC);
    FlinwaveHalfRight(i)=trapz(x,dphidtS.*nu3CR);
    FnonlinwaveHalfLeft(i)=trapz(x,deletaKS.*nu3LC);
    FnonlinwaveHalfRight(i)=trapz(x,deletaKS.*nu3CR);
    nu1LC        = nux;nu1CR        = nux;
    nu1LC(IndxC+1:IndxR)=0;nu1CR(IndxL:IndxC)=0;
    FlinwaveHalfLeft(Nship+i)=trapz(x,dphidtS.*nu1LC);
    FlinwaveHalfRight(Nship+i)=trapz(x,dphidtS.*nu1CR);
    FnonlinwaveHalfLeft(Nship+i)=trapz(x,deletaKS.*nu1LC);
    FnonlinwaveHalfRight(Nship+i)=trapz(x,deletaKS.*nu1CR);
end
end