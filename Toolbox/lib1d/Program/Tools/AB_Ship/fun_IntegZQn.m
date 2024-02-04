function IZDQn = fun_IntegZQn(kn,eta,D,z1,z2)
y2=D+z2;
y1=D+z1;
IZDQn=(z2.*sin(kn.*y2)./kn+cos(kn.*y2)/kn/kn-(z1.*sin(kn.*y1)./kn+cos(kn.*y1)/kn/kn))./cos(kn.*(D+eta));