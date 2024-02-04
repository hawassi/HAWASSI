function IZD2varphin = fun_IntegZD2Varphin(mun,D,z1,z2)
y2=D+z2;
y1=D+z1;
if mun~=0
IZD2varphin=(sin(mun.*y2).*(y2^2-2/mun/mun)+2*y2.*cos(mun*y2)/mun-(sin(mun.*y1).*(y1^2-2/mun/mun)+2*y1.*cos(mun*y1)/mun))./mun;
else
IZD2varphin=(y2^3-y1^3)/3;    
end