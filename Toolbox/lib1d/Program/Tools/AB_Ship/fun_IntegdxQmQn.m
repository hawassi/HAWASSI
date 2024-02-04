function Int1mn=fun_IntegdxQmQn(m,n,km,kn,eta,D,dxeta,dxD,Z)
%integral_-D_Zmin Q_m Qndz
H=D+eta;
dxH=dxD+dxeta;
if m==n
    Int1mn=km.*(dxD.*(cos(2*km*(D+Z))-1)/4/km+dxH.*tan(km.*H)/2.*((D+Z)+sin(2*km*(D+Z))/2/km))./cos(km.*H).^2;
else
    A=(km+kn);
    B=(km-kn); 
    Int1mn=(-dxD.*km.*(-(cos(A.*(D+Z))/A+cos(B.*(D+Z))/B)+2)/2+dxH.*km.*tan(km.*H).*(sin(A.*(D+Z))/A+sin(B.*(D+Z))/B)/2)/cos(km*H)/cos(kn*H);
end
end