function Cp2Var_hat=OpCp2_var_Wall_B(Rho_k,var_hat,indWi,indWf,Oprt)
var        =Ifft(var_hat);
if Oprt.interp==2
Cp2Var_hat  = HSS_Wall_k(Oprt.Csqmin,Oprt.gammin,var_hat,var,Rho_k,indWi,indWf)+...
             HSS_Wall_k(Oprt.Csqplus,Oprt.gamplus,var_hat,var,Rho_k,indWi,indWf);
else
Cp2Var_hat  = HSS_Wall_k(Oprt.Csqmin,Oprt.gammin,var_hat,var,Rho_k,indWi,indWf)+...
              HSS_Wall_k(Oprt.Csqmid,Oprt.gammid,var_hat,var,Rho_k,indWi,indWf)+...
             HSS_Wall_k(Oprt.Csqplus,Oprt.gamplus,var_hat,var,Rho_k,indWi,indWf);    
end
