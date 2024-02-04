function M1phi_hat=funOp_M1(L0,eta,k,phi)
M1phi_hat=-L0.*fft(eta.*ifft(L0.*fft(phi),'symmetric'))-...
    1i.*k.*fft(eta.*ifft(1i.*k.*fft(phi),'symmetric'));
end