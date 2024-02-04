function div=funOprt_div2d(Kx,Ky,func_hat)
div=1i.*Kx.*func_hat.x_hat+1i.*Ky.*func_hat.y_hat;
