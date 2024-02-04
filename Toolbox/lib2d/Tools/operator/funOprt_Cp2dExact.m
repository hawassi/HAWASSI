function Cp2d=funOprt_Cp2dExact(Kx,Ky,depth,g)
K=sqrt(Kx.^2+Ky.^2);
Cp2d=sqrt(g*tanh(K.*depth)./K);
Cp2d(K==0)=sqrt(g.*depth);