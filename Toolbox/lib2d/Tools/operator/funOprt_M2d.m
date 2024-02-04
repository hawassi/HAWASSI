function M2dA_hat=funOprt_M2d(g,Kx,Ky,Cp2dSqInterp,u_hat,u,v_hat,v)

Cp2gradphi_hat.x_hat = funOprt_Cp2d(Cp2dSqInterp,u_hat,u);
Cp2gradphi_hat.y_hat = funOprt_Cp2d(Cp2dSqInterp,v_hat,v);
divCp2gradphi_hat    = funOprt_div2d(Kx,Ky,Cp2gradphi_hat);
M2dA_hat             = -divCp2gradphi_hat/g;   
