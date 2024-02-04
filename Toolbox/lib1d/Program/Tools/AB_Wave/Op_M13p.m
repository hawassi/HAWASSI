function M1 = Op_M13p(k,g,Cm,gm,Cp,gp,Cc,gc,A)
iftA=Ifft(A);
   HSSs_hat=HSS(Cm,gm,A,iftA)+HSS(Cp,gp,A,iftA)+HSS(Cc,gc,A,iftA);
   M1=-(HSSs_hat).*1i.*k./g; 
end