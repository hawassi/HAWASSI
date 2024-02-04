function C2u_hat=OpC2interp(Oprt,u_hat)
u=Ifft(u_hat);
if Oprt.interp==2
C2u_hat=HSS(Oprt.Csqmin,Oprt.gammin,u_hat,u)+...
        HSS(Oprt.Csqplus,Oprt.gamplus,u_hat,u);
else
C2u_hat=HSS(Oprt.Csqmin,Oprt.gammin,u_hat,u)+...
        HSS(Oprt.Csqmid,Oprt.gammid,u_hat,u)+...
        HSS(Oprt.Csqplus,Oprt.gamplus,u_hat,u);
end
end