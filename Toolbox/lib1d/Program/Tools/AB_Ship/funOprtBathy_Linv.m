function Linvphi_hat=funOprtBathy_Linv(g,kk,Oprt,eta,phi,OrdeOpL)
L0invphi_hat=funOprtBathy_L0inv(g,kk,Oprt,phi);
Linvphi_hat=L0invphi_hat;
if OrdeOpL~=1
dxetadxL0invphi=Ifft(1i.*kk.*fft(eta.*Ifft(1i.*kk.*L0invphi_hat)));
L0invdxetadxL0invphi_hat=funOprtBathy_L0inv(g,kk,Oprt,dxetadxL0invphi);  
Linvphi_hat=Linvphi_hat+fft(eta.*phi)+L0invdxetadxL0invphi_hat;
if OrdeOpL~=2
dxxeta2phi=Ifft((1i.*kk).^2.*fft(eta.^2.*phi));
L0invdxxeta2phi_hat=funOprtBathy_L0inv(g,kk,Oprt,dxxeta2phi); 
dxetadxetaphi=Ifft(1i.*kk.*fft(eta.*ifft(1i.*kk.*fft(eta.*phi))));
L0invdxetadxetaphi_hat=funOprtBathy_L0inv(g,kk,Oprt,dxetadxetaphi); 
dxetadxLinvphi=Ifft(1i.*kk.*fft(eta.*Ifft(1i.*kk.*L0invphi_hat)));
L0invdxetadxLinvphi_hat=funOprtBathy_L0inv(g,kk,Oprt,dxetadxLinvphi); 
dxetadxL0invdxetadxLinvphi=Ifft(1i.*kk.*fft(eta.*Ifft(1i.*kk.*L0invdxetadxLinvphi_hat)));
L0invdxetadxL0invdxetadxLinvphi_hat=funOprtBathy_L0inv(g,kk,Oprt,dxetadxL0invdxetadxLinvphi); 

term3order=-0.5.*(L0invdxxeta2phi_hat+fft(eta.^2.*Ifft((1i.*kk).^2.*L0invphi_hat)))...
    +fft(eta.*Ifft(1i.*kk.*fft(eta.*Ifft(1i.*kk.*L0invphi_hat))))...
    +L0invdxetadxetaphi_hat+...
    +L0invdxetadxL0invdxetadxLinvphi_hat;
Linvphi_hat=Linvphi_hat+term3order;
end
end
end