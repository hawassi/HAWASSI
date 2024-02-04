function Fphi_hat=funOprt_Flux(g,kk,Cp,eta,phi,OrdeOp)
F0=Cp.^2.*1i.*kk/g;
phi_hat=fft(phi);
Fphi_hat=F0.*phi_hat;
if OrdeOp~=1
    F1phi_hat=funOp_F1(F0,eta,kk,phi);
    Fphi_hat=Fphi_hat+ F1phi_hat;
    if OrdeOp~=2
        F2phi_hat=funOp_F2(F0,eta,kk,phi);
        Fphi_hat=Fphi_hat+ F2phi_hat;
    end
end