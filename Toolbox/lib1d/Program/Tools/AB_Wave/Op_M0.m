function M0A_hat = Op_M0(k,g,Oprt,A_hat)
A=Ifft(A_hat);
if Oprt.interp==2
M0A_hat=-1i.*k.*(HSS(Oprt.Csqmin,Oprt.gammin,A_hat,A)+...
    HSS(Oprt.Csqplus,Oprt.gamplus,A_hat,A))/g; 
else
M0A_hat=-1i.*k.*(HSS(Oprt.Csqmin,Oprt.gammin,A_hat,A)+...
    HSS(Oprt.Csqmid,Oprt.gammid,A_hat,A)+...
    HSS(Oprt.Csqplus,Oprt.gamplus,A_hat,A))/g;      
end
end