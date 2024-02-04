function Heav2d=funC_Heav2d(Y,b,t,Nx)


if b==0
Heavb=1;
else 
Heavb =funC_Heaviside(Y-(Y(1)+b));
end

if t==0
Heavt=1;
else 
Heavt = (1-funC_Heaviside(Y-(Y(end)-t)));
end

Heav=Heavb.*Heavt;   

if b==0 && t==0
Heav2d   =1; 
else
Heav2d   = repmat(Heav',1,Nx);
end
