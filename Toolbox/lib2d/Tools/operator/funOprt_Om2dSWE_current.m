function om2d = funOprt_Om2dSWE_current(Kx,Ky,depth,g,~,ux,uy) 
K=sqrt(Kx.^2+Ky.^2);
om2d= sqrt(g*depth).*K+Kx.*ux+Ky.*uy; 