function L0invphi_hat=funOprtBathy_L0inv(g,kk,Oprt,phi)
L0min     =-1i.*kk.*Oprt.Csqmin.*1i.*kk./g;
L0plus     =-1i.*kk.*Oprt.Csqplus.*1i.*kk./g;
L0invmin  = 1./L0min;
L0invplus  = 1./L0plus;
L0invmin(abs(L0invmin)==Inf)=0;
L0invplus(abs(L0invplus)==Inf)=0;
gm          = Oprt.gammin;
gp          = Oprt.gamplus;

L0invphi_hat   = L0invplus.*fft(gp.*phi)+L0invmin.*fft(gm.*phi);%HSS(L0invplus,gp,phi_hat,phi);
L0invphi_hat(1)=0;
