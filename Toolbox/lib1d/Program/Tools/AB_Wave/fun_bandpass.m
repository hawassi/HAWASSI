function sig=fun_bandpass(timesig,insig,LFreq,HFreq)
sighat=zeros(size(insig));
ww=funC_freqspace(timesig);
insighat=fft(insig);
indW1p1=closest(ww,LFreq);
indW1p2=closest(ww,HFreq);

indW1m1=closest(ww,-HFreq);
indW1m2=closest(ww,-LFreq);
sighat(indW1p1:indW1p2)=insighat(indW1p1:indW1p2);
sighat(indW1m1:indW1m2)=insighat(indW1m1:indW1m2);
sig=Ifft(sighat);
