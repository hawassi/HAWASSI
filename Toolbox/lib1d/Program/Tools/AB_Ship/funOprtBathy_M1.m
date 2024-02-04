function [M1phi_hat,LetaLphi_hat]=funOprtBathy_M1(g,Oprt,eta,kk,phi,Lphi_hat)
etaLphi_hat=fft(eta.*Ifft(Lphi_hat));
etaLphi=Ifft(etaLphi_hat);
LetaLphi_hat=funOprtBathy_L0(g,kk,Oprt,etaLphi);
M1phi_hat=-LetaLphi_hat-...
    1i.*kk.*fft(eta.*Ifft(1i.*kk.*fft(phi)));
end