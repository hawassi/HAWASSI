function M2phi_hat=funOprtBathy_M2(g,Oprt,eta,k,phi,Lphi_hat,LetaLphi_hat)

Lphi=ifft(Lphi_hat,'symmetric');
LetaLetaLphi_hat=funOprtBathy_L0(g,k,Oprt,eta.*Ifft(LetaLphi_hat));
M2phi_hat=LetaLetaLphi_hat+...
    0.5*(1i.*k.*1i.*k.*fft(eta.^2.*Lphi)+...
    funOprtBathy_L0(g,k,Oprt,eta.^2.*Ifft(1i.*k.*1i.*k.*fft(phi))));
end