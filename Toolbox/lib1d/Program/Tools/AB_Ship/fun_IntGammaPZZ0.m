function IntGammaP=fun_IntGammaPZZ0(D,Zc,X,Z0,Z1)
Z=[Z1 Z0];
IntGammaPvect=Z.^4./4+2*Z.^3*D/3+D^2*Z.^2/2-(Z+D).^3*Zc*X/3-X^3*Z.^2/6+X^3*Zc.*Z./3;
IntGammaP=(IntGammaPvect(1)-IntGammaPvect(2))/2/(D+Z1);