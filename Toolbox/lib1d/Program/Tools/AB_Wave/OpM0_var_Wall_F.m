function M0var_hat = OpM0_var_Wall_F(k,g,Oprt,Rho_k,var_hat,indWi,indWf)
Cp2      =Oprt.Upd.^2;
RhoCp2Var =Ifft(Cp2.*Rho_k.*var_hat);
Cp2Var    =Ifft(Cp2.*var_hat);
Cp2Var(indWi:indWf)=RhoCp2Var(indWi:indWf);
Cp2Var_hat=fft(Cp2Var);
M0var_hat  =-1/g.*1i.*k.*Cp2Var_hat; 