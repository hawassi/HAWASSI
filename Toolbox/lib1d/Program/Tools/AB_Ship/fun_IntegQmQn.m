function Int1mn=fun_IntegQmQn(m,n,kappa_m,kappa_n,eta,D,Z)
%integral_-D_Zmin Q_m Qndz
if m==n
    Int1mn=(D+Z+sin(2*kappa_m*(D+Z))/2/kappa_m)/(2*cos(kappa_m*(D+eta))^2);
else
    A=(kappa_m+kappa_n);
    B=(kappa_m-kappa_n); 
    Int1mn=(sin(A*(D+Z))/A+sin(B*(D+Z))/B)/2/cos(kappa_m*(D+eta))/cos(kappa_n*(D+eta));
end
end