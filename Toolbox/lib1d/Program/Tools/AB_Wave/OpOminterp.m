function OmVar_hat=OpOminterp(Oprt,var_hat)
var=Ifft(var_hat);
if Oprt.interp==2
OmVar_hat=HSS(Oprt.Ommin,Oprt.gammin,var_hat,var)+...
        HSS(Oprt.Omplus,Oprt.gamplus,var_hat,var);
else
OmVar_hat=HSS(Oprt.Ommin,Oprt.gammin,var_hat,var)+...
        HSS(Oprt.Ommid,Oprt.gammid,var_hat,var)+...
        HSS(Oprt.Omplus,Oprt.gamplus,var_hat,var);
end
end