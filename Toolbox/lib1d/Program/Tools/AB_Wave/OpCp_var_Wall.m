function CpVar_hat=OpCp_var_Wall(Cp,Rho_k,var_hat,indWi,indWf)
RhoCpVar=Ifft(Cp.*Rho_k.*var_hat);
CpVar   =Ifft(Cp.*var_hat);
CpVar(indWi:indWf)=RhoCpVar(indWi:indWf);
CpVar_hat=fft(CpVar);