function [Lphi_hat,FluxPhi_hat]=funOprt_Lbathy(g,k,Oprt,eta,u,OrdeOpL,ChiAdj)
u_hat=fft(u);
FluxPhi_hat=OpC2interp(Oprt,u_hat)./g;
Lphi_hat= -1i*k.*FluxPhi_hat;
if OrdeOpL~=1
    M0u_hat      = Op_M0(k,g,Oprt,u_hat);
    M0u          = Ifft(M0u_hat);
    etaM0u_hat   = fft(eta.*M0u);
    M1etaM0u_hat = Op_M1(k,g,Oprt,etaM0u_hat);
    FluxPhi2_hat  =fft(eta.*u)+M1etaM0u_hat;
    Lphi_hat     = Lphi_hat-fft(ChiAdj.*Ifft(1i.*k.*(FluxPhi2_hat)));
    FluxPhi_hat=FluxPhi_hat+FluxPhi2_hat;
    if OrdeOpL~=2
        disp('third order is not available yet.')
    end
end