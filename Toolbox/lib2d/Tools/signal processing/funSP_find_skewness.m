function Sk=funSP_find_skewness(etaiter)
meaneta=mean(etaiter);
Sk=mean((etaiter-meaneta).^3)/((mean((etaiter-meaneta).^2))^(3/2));    