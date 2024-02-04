function Lphi_hat=funOprtBathy_L(g,kk,Oprt,eta,phi,OrdeOpL,ChiAdj)

Lphi_hat=funOprtBathy_L0(g,kk,Oprt,phi);
if OrdeOpL~=1
L0phi_hat=Lphi_hat;
[M1phi_hat,L0etaL0phi_hat]=funOprtBathy_M1(g,Oprt,eta,kk,phi,L0phi_hat);
Lphi_hat=Lphi_hat+ fft(Ifft(M1phi_hat).*ChiAdj);  
if OrdeOpL~=2
M2phi_hat=funOprtBathy_M2(g,Oprt,eta,kk,phi,L0phi_hat,L0etaL0phi_hat);
Lphi_hat=Lphi_hat+ fft(Ifft(M2phi_hat).*ChiAdj);     
end
end