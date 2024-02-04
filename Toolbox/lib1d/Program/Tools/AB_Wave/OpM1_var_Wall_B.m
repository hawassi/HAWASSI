function M1Var_hat=OpM1_var_Wall_B(k,g,Rho_k,var_hat,indWi,indWf,Oprt)
dxvar_hat    =1i.*k.*var_hat;
Cp2dxVar_hat=OpCp2_var_Wall_B(Rho_k,dxvar_hat,indWi,indWf,Oprt);
M1Var_hat   = -Cp2dxVar_hat/g;  