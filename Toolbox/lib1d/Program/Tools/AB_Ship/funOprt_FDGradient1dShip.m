function dxshape=funOprt_FDGradient1dShip(shape,x,Sxc)
dx=x(2)-x(1);
dxshapeF=funOprt_FDGradient1d(shape,dx,1);
dxshapeB=funOprt_FDGradient1d(shape,dx,-1);
dxshapeB(x<Sxc)=0;
dxshapeF(x>=Sxc)=0;
dxshape=(dxshapeF+dxshapeB);