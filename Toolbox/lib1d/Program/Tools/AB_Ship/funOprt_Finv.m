function Finv_varhat=funOprt_Finv(k,F0,eta,varhat,OrdeOpF)
L0=-1i.*k.*F0;
F0inv=1./F0;
Finv_varhat=F0inv.*varhat;
if OrdeOpF~=1
F0inv_varhat= F0inv.*varhat;   
L0F0inv_varhat=L0.*F0inv_varhat;
L0F0inv_varhat(isnan(L0F0inv_varhat))=0;
dxF0inv_varhat=1i.*k.*F0inv_varhat;
dxF0inv_varhat(isnan(dxF0inv_varhat))=0;
F1F0inv_varhat=F0.*fft(eta.*Ifft(L0F0inv_varhat))+fft(eta.*Ifft(dxF0inv_varhat)); 
F0invF1F0inv_varhat=F0inv.*F1F0inv_varhat;
Finv_varhat=Finv_varhat-F0invF1F0inv_varhat;%+(F0invF1F0inv_varhat).^2;    
end
Finv_varhat(1)=0;
