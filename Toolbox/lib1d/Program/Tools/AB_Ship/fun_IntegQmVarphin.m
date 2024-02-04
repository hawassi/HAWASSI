function Int3mn=fun_IntegQmVarphin(m,n,kappa_m,mu_n,eta,D,Z)
%integral_z1_z2 Q_m varphi_n dz
if m==0 && n==0
  Int3mn=sin(kappa_m*(D+Z))/(kappa_m*cos(kappa_m*(D+eta)));
else
  Int3mn=kappa_m*sin(kappa_m*(D+Z))*cos(mu_n*(D+Z))/cos(kappa_m*(D+eta))/(kappa_m^2-mu_n^2);  
end
end