function IZD2Qn = fun_IntegZD2Qn(kn,eta,D,z1,z2)
y2=D+z2;
y1=D+z1;
IZD2Qn=(sin(kn.*y2).*(y2^2-2/kn/kn)+2*y2.*cos(kn*y2)/kn-(sin(kn.*y1).*(y1^2-2/kn/kn)+2*y1.*cos(kn*y1)/kn))./kn./cos(kn.*(D+eta));