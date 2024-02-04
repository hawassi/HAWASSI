function [Fev,Lev]=fun_wall_condition(k,u,eta,parwall)
Fev=zeros(size(u));
Lev=zeros(size(u));
indWl=parwall.indx(1);
indWr=parwall.indx(2);

phi=Ifft(fft(u)./1i./k);
phiL=phi(indWl);
phiR=phi(indWr);
uL=u(indWl);%
uR=u(indWr);%

[AmatW,BvectW] =fun_construct_matrices_reflecting_wall(parwall.Dwall,...
       parwall.Evmodes,parwall.kappa,parwall.alphaR,phiL,uL);

CoefW=AmatW\BvectW