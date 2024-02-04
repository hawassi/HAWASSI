function L0phi_hat=funOprtBathy_L0(g,kk,Oprt,phi)
Cm2         = Oprt.Csqmin;
Cp2         = Oprt.Csqplus;
gm          = Oprt.gammin;
gp          = Oprt.gamplus;
u_hat= 1i.*kk.*fft(phi);u=Ifft(u_hat);
Csq_u_hat   = HSS(Cm2,gm,u_hat,u)+HSS(Cp2,gp,u_hat,u);
L0phi_hat   = -1i*kk.*Csq_u_hat/g;