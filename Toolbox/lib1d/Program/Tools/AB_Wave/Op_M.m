function M = Op_M(k,g,Cm,gm,Cp,gp,A)
iftA=Ifft(A);
    M=-1i.*k.*(HSS(Cm,gm,A,iftA)+HSS(Cp,gp,A,iftA))/g; 
end