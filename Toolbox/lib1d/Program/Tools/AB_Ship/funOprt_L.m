function Lphi_hat=funOprt_L(kk,L0,eta,phi,OrdeOpL,ChiAdj)
phi_hat=fft(phi);
Lphi_hat=L0.*phi_hat;
if OrdeOpL~=1
M1phi_hat=funOp_M1(L0,eta,kk,phi);
Lphi_hat=Lphi_hat+ fft(Ifft(M1phi_hat).*ChiAdj);  
if OrdeOpL~=2
M2phi_hat=funOp_M2(L0,eta,kk,phi);
Lphi_hat=Lphi_hat+ fft(Ifft(M2phi_hat).*ChiAdj);     
end
end

