function tauU=fun_tangential_potential_anzats(k,x,L0,F0,taug,eta,phi,Nonlin,ChiAdj,xShip,chiSWl,Ns)

% if Nonlin<=2
%    Nonlin=2; 
% end
mu=-Ifft(funOprt_L(k,L0,eta,phi,Nonlin,ChiAdj));
mubar=mu;
dx=x(2)-x(1);
for i=1:Ns
    SL=xShip(i,3)-xShip(i,1);
   indxL=closest(x,xShip(i,1)-5*SL);
   indxR=closest(x,xShip(i,3)+5*SL);
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

phidifbar=Ifft(funOprt_Linv(k,L0,eta,mubar,Nonlin));
phidif=zeros(size(phidifbar));
for i=1:Ns
indxL=closest(x,xShip(i,1))-1; % due to the barge make sure the flux at ship bottom at water line is not zero.
indxR=closest(x,xShip(i,3))+1;
A=(phidifbar(indxR).*taug(indxL)-phidifbar(indxL).*taug(indxR))/(taug(indxR)-taug(indxL));
B=(phidifbar(indxR)-phidifbar(indxL))/(taug(indxL)-taug(indxR));
phidif(indxL:indxR)=phidifbar(indxL:indxR)+A+B.*taug(indxL:indxR);
end
tauU=phi+phidif;
% figure;
% subplot(3,1,1);
% plot(x,mu,'r',x,mubar,'b')
% xlabel('x [m]');legend('\mu','\bar{\mu}')
% subplot(3,1,2);
% plot(x,phi,'r',x,phidif,'b',x,tauU,'k',x,phidifbar,'g')
% xlabel('x [m]');legend('\phi','\phi_d','\tau')
% subplot(3,1,3)
% Lphidif=Ifft(funOprt_L(k,L0,eta,phidif,Nonlin,ChiAdj));
% Ltau0=Ifft(funOprt_L(k,L0,eta,tauU,Nonlin,ChiAdj));
% plot(x,-mu,'r',x,Lphidif,'b',x,Ltau0,'k')
% xlabel('x [m]');legend('L\phi','L\phi_d','L\tau')