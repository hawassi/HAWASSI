function Int1mn=fun_IntegdxVarphimQn(mu_m,kn,eta,D,dxD,Z_min)
%integral_-D_Zmin dxVarphi_m Qndz
H=D+eta;
A=(mu_m+kn);
B=(mu_m-kn); 
Int1mn=-dxD.*mu_m.*(cos(A.*(D+Z_min))./A+cos(B.*(D+Z_min))./B-(1./A+1./B))./2./cos(kn.*H);
end