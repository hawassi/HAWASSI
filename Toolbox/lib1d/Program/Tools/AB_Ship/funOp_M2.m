function M2phi_hat=funOp_M2(L0,eta,k,phi)
Lphi_hat=L0.*fft(phi);
Lphi=ifft(Lphi_hat,'symmetric');
M2phi_hat=L0.*fft(eta.*ifft(L0.*fft(eta.*Lphi),'symmetric'))+...
    0.5*(1i.*k.*1i.*k.*fft(eta.^2.*Lphi)+L0.*fft(eta.^2.*ifft(1i.*k.*1i.*k.*fft(phi),'symmetric')));
end