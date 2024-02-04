% HybridSpatialSpectral operator; 
% from fourier function z to fourier function Az; a is space interp function
% iftz = ifft(z,'symmetric'); requires 3 FT's

function Az_hat = HSS_Wall_k(A,a,z_hat,z,Rho_k,indWi,indWf)

RhoAz          =Ifft(A.*Rho_k.*z_hat);
Az             =Ifft(A.*z_hat);
Az(indWi:indWf)=RhoAz(indWi:indWf);
Az_hat         =fft(Az);

az_hat          =fft(a.*z); 
RhoAaz          =Ifft(A.*Rho_k.*az_hat);
Aaz(indWi:indWf)=RhoAaz(indWi:indWf);      
Aaz_hat         =fft(Aaz)';
Az_hat = ( fft(a.*Ifft(Az_hat)) + Aaz_hat)/2;  
end