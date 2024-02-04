function Fphi_hat=funOprt_FluxWall(g,kk,Cp,eta,phi,OrdeOp,Wallchar)
phi_hat=fft(phi);
u_hat=1i.*kk.*phi_hat;
Cu_hat          = Cp.*u_hat;
Cu              = Ifft(Cu_hat);
Cu_Wall         = Cu.*Wallchar; 

Fphi_hat=Cp.*fft(Cu_Wall)./g;


% if OrdeOp~=1
%     F1phi_hat=funOp_F1(F0,eta,kk,phi);
%     Fphi_hat=Fphi_hat+ F1phi_hat;
%     if OrdeOp~=2
%         F2phi_hat=funOp_F2(F0,eta,kk,phi);
%         Fphi_hat=Fphi_hat+ F2phi_hat;
%     end
% end