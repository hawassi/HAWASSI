function F1phi_hat=funOp_F1(F0,eta,k,phi)
F1phi_hat=-F0.*fft(eta.*ifft(-1i*k.*F0.*fft(phi),'symmetric'))+...
        fft(eta.*ifft(1i.*k.*fft(phi),'symmetric'));
end