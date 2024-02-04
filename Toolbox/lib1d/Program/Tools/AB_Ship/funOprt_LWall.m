function Lphi_hat=funOprt_LWall(g,kk,Cp,eta,phi,OrdeOpL,ChiAdj,Wallchar)
phi_hat=fft(phi);
u_hat=1i.*kk.*phi_hat;
Cu_hat          = Cp.*u_hat;
Cu              = Ifft(Cu_hat);
Cu_Wall         = Cu.*Wallchar; 

Lphi_hat=-1i.*kk.*Cp.*fft(Cu_Wall)./g;
% if OrdeOpL~=1
% M1phi_hat=funOp_M1(L0,eta,kk,phi);
% Lphi_hat=Lphi_hat+ fft(Ifft(M1phi_hat).*ChiAdj);  
% if OrdeOpL~=2
% M2phi_hat=funOp_M2(L0,eta,kk,phi);
% Lphi_hat=Lphi_hat+ fft(Ifft(M2phi_hat).*ChiAdj);     
% end
end