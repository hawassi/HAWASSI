function Int3mn=fun_diff_Integ3mn(m,n,kappa_m,mu_n,D,T)
%integral_z1_z2 chi_m chi_n dz
if m==0 && n==0
    Int3mn=sin(kappa_m*(D-T))/(kappa_m*cos(kappa_m*D));
else
  Int3mn=(-1)^(n)*kappa_m*sin(kappa_m*(D-T))/cos(kappa_m*D)/(kappa_m^2-mu_n^2);  
end
end