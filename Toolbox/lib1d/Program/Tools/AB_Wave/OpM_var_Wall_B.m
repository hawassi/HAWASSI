function MVar_hat=OpM_var_Wall_B(k,g,Rho_k,var_hat,indWi,indWf,Oprt)
Cm2        = Oprt.Csqmin;
Cp2        = Oprt.Csqplus;
gm         = Oprt.gammin;
gp         = Oprt.gamplus;
var        =Ifft(var_hat);
Cp2Var_hat = HSS_Wall_k(Cm2,gm,var_hat,var,Rho_k,indWi,indWf)+...
             HSS_Wall_k(Cp2,gp,var_hat,var,Rho_k,indWi,indWf);
MVar_hat   = -1i.*k.*Cp2Var_hat/g;              