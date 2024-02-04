function [Hs,Tp,Tm01]=funSP_find_stat_Hs_Tp_Tm01_fromSignal(time,Eta,Spsmooth)

ww=funC_freqspace(time);
Eta=Eta-mean(Eta);
varSig=var(Eta);
S_hat=abs(fft(Eta)).^2;
var0=trapz(ww(1:floor(end/2)),S_hat(1:floor(end/2)));
OmS=ww(1:floor(end/2));
S_hat=S_hat(1:floor(end/2));
S_hat=S_hat*varSig/var0;  % variance density spectrum

if Spsmooth>1
     S_hat=smooth(S_hat,Spsmooth) ;  
end

[~,sppeakindexPP]=max(S_hat);
nupeak     =ww(sppeakindexPP);   
Hs      =4.*sqrt(varSig);
Tm01    =2*pi*trapz(OmS,S_hat)/trapz(OmS,OmS.*S_hat);
Tp      = 2*pi/nupeak;