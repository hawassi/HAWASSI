function Ug2d=funOprt_Ug2d(Kx,Ky,depth,g,ej)
K=sqrt(Kx.^2+Ky.^2);
Ug2d = sign(Kx*ej(1)+Ky*ej(2)).*sqrt(g)/2./(K.*tanh(depth.*K)).^(1/2).*(tanh(depth*K)+K.*(1-tanh(depth*K).^2)*depth);
Ug2d(K==0)=sqrt(g*depth);