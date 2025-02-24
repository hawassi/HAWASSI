function [OmS,S_hat]=funSP_variance_density_spectrum1d(time,Eta)
omsig=funC_freqspace(time);
Eta=Eta-mean(Eta);
varSig=var(Eta);
S_hat=abs(fft(Eta)).^2;
var0=trapz(omsig(1:floor(end/2)),S_hat(1:floor(end/2)));
OmS=omsig(1:floor(end/2));
S_hat=S_hat(1:floor(end/2));
S_hat=S_hat*varSig/var0;