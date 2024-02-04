function Om2d=funOprt_Om2dExact(Kx,Ky,depth,g,ej)
K=sqrt(Kx.^2+Ky.^2);
Om2d=sign(Kx*ej(1)+Ky*ej(2)).*sqrt(g*K.*tanh(K.*depth));
