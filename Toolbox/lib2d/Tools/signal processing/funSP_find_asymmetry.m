function As=funSP_find_asymmetry(etaiter)
meaneta=mean(etaiter);
etaNow=etaiter-meaneta;
Hilbeta=imag(hilbert(etaNow));
As=mean(Hilbeta.^3)/((mean((etaiter-meaneta).^2))^(3/2));