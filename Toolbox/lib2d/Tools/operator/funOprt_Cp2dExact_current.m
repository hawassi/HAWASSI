function Cp2d=funOprt_Cp2dExact_current(Kx,Ky,depth,g,ux,uy)
K=sqrt(Kx.^2+Ky.^2);
Cp2d=sqrt(g*tanh(K.*depth)./K)+(Kx.*ux+Ky.*uy)./K;
Cp2d(K==0)=sqrt(g.*depth);