function F2phi_hat=funOp_F2(F0,eta,k,phi)
Lphi_hat=-1i.*k.*F0.*fft(phi);
Lphi=ifft(Lphi_hat,'symmetric');
F2phi_hat=F0.*fft(eta.*ifft(-1i.*k.*F0.*fft(eta.*Lphi),'symmetric'))+...
    0.5*(-1i.*k.*fft(eta.^2.*Lphi)+F0.*fft(eta.^2.*ifft(1i.*k.*1i.*k.*fft(phi),'symmetric')));
end