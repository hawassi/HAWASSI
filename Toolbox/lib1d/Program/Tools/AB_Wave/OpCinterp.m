function Cu_hat=OpCinterp(Oprt,u_hat)
u=Ifft(u_hat);
if Oprt.interp==2
Cu_hat=HSS(sqrt(Oprt.Csqmin),Oprt.gammin,u_hat,u)+...
       HSS(sqrt(Oprt.Csqplus),Oprt.gamplus,u_hat,u);
else
Cu_hat=HSS(sqrt(Oprt.Csqmin),Oprt.gammin,u_hat,u)+...
       HSS(sqrt(Oprt.Csqmid),Oprt.gammid,u_hat,u)+...
       HSS(sqrt(Oprt.Csqplus),Oprt.gamplus,u_hat,u);   
end
end