function M0 = Op_M03p(k,g,Cm,gm,Cp,gp,Cc,gc,A)
iftA=Ifft(A);
    HSSs_hat=HSS(Cm,gm,A,iftA)+HSS(Cp,gp,A,iftA)+HSS(Cc,gc,A,iftA);
    M0=-1i.*k.*(HSSs_hat)/g; 
end