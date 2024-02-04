function IVarphin = fun_IntegVarphin(mun,D,z1,z2)
if mun~=0
IVarphin=(sin(mun.*(D+z2))-sin(mun.*(D+z1)))/mun;
else
IVarphin=z2-z1;    
end