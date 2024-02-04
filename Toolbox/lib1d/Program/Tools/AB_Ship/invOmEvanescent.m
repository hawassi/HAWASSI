function kappah = invOmEvanescent(sig0D,modes)

kappah=zeros(modes,1);
iter=0;
for ii=1:modes    
kappah(ii)= fzero(@(kappah) -sig0D/kappah-tan(kappah),[(ii+iter)*pi/2+0.01 (ii+iter+2)*pi/2-0.01]); 
iter=iter+1;
end