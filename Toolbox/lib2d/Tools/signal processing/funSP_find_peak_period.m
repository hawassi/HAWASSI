function Tp=funSP_find_peak_period(T,etaiter)
ww=funC_freqspace(T);
omsig=ww(1:floor(end/2));
[~,sppeakindex]=max(abs(fft(etaiter)));
 nupeak      = abs(omsig(sppeakindex)); % this is peak-frequency!!!!
 Tp=2*pi/nupeak;      