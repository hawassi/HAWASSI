function M0Var_hat=OpM0_var_Wall_B(k,g,Rho_k,var_hat,indWi,indWf,Oprt)
Cp2Var_hat=OpCp2_var_Wall_B(Rho_k,var_hat,indWi,indWf,Oprt);
M0Var_hat   = -1i.*k.*Cp2Var_hat/g;              