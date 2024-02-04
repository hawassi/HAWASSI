function Linvphi_hat=funOprt_Linv(kk,L0,eta,phi,OrdeOpL)
L0inv=1./L0;
phihat=fft(phi);
Linvphi_hat=L0inv.*phihat;



if OrdeOpL~=1
term2order=fft(eta.*phi)+L0inv.*1i.*kk.*fft(eta.*Ifft(1i.*kk.*L0inv.*phihat));
Linvphi_hat=Linvphi_hat+term2order;
if OrdeOpL~=2
term3order=-0.5.*(L0inv.*(1i.*kk).^2.*fft(eta.^2.*phi)+fft(eta.^2.*Ifft((1i.*kk).^2.*L0inv.*phihat)))...
    +fft(eta.*Ifft(1i.*kk.*fft(eta.*Ifft(1i.*kk.*L0inv.*phihat))))...
    +L0inv.*1i.*kk.*fft(eta.*Ifft(1i.*kk.*fft(eta.*phi)))+...
    +L0inv.*1i.*kk.*fft(eta.*Ifft(1i.*kk.*L0inv.*1i.*kk.*fft(eta.*Ifft(1i.*kk.*L0inv.*phihat))));
Linvphi_hat=Linvphi_hat+term3order;
end
end
Linvphi_hat(1)=0;
end