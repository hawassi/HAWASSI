function IntdxQm=fun_IntegdxQn(k,eta,D,dxeta,dxD,z1,z2)
H=D+eta;
dxH=dxD+dxeta;
IntdxQm=(sin(k.*(D+z2))-sin(k.*(D+z1))).*sin(k.*H).*dxH./cos(k.*H).^2+(cos(k.*(D+z2))-cos(k.*(D+z1))).*dxD.*k./cos(k.*H);

