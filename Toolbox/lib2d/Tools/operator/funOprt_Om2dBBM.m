function Om2d = funOprt_Om2dBBM(Kx,Ky,depth,g,~) 
K=sqrt(Kx.^2+Ky.^2);
Om2d=sqrt(g.*depth)*K./(1+(depth*K).^2/6);

