function IZVarphin = fun_IntegZVarphin(mu_n,D,z1,z2)
if mu_n==0
IZVarphin=z2^2/2-(z1^2/2);    
else
IZVarphin=z2*sin(mu_n*(D+z2))/mu_n+cos(mu_n*(D+z2))/mu_n/mu_n ...
           -(z1*sin(mu_n*(D+z1))/mu_n+cos(mu_n*(D+z1))/mu_n/mu_n) ;
end