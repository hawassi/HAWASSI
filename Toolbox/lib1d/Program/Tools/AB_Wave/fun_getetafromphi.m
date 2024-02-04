function eta=fun_getetafromphi(g,phi,Om) 
etahat=fft(phi).*(1i.*Om)/g;
eta=Ifft(etahat);