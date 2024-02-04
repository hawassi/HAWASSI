function om2d = funOprt_Om2dKdV(Kx,Ky,depth,g,~) 
K=sqrt(Kx.^2+Ky.^2);
om2d= sqrt(g*depth)*K.*(1-(depth*K).^2/6); 

