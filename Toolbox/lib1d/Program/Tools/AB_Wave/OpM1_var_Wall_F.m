function M1var_hat = OpM1_var_Wall_F(k,g,Oprt,Rho_k,var_hat,indWi,indWf)
Cp2         =Oprt.Upd.^2;
dxvar_hat   =1i.*k.*var_hat;
RhoCp2dxVar =Ifft(Cp2.*Rho_k.*dxvar_hat);
Cp2dxVar    =Ifft(Cp2.*dxvar_hat);
Cp2dxVar(indWi:indWf)=RhoCp2dxVar(indWi:indWf);
Cp2dxVar_hat=fft(Cp2dxVar);
M1var_hat   =-1/g.*Cp2dxVar_hat; 