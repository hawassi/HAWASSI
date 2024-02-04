function M1A_hat = Op_M1(k,g,Oprt,A_hat)
A=Ifft(A_hat);
if Oprt.interp==2
M1A_hat=-(HSS(Oprt.Csqmin,Oprt.gammin,A_hat,A)+...
    HSS(Oprt.Csqplus,Oprt.gamplus,A_hat,A)).*1i.*k/g; 
else
M1A_hat=-(HSS(Oprt.Csqmin,Oprt.gammin,A_hat,A)+...
    HSS(Oprt.Csqmid,Oprt.gammid,A_hat,A)+...
    HSS(Oprt.Csqplus,Oprt.gamplus,A_hat,A)).*1i.*k/g;      
end

end