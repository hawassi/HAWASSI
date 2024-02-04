function grad=funOprt_grad2d_hat(Kx,Ky,func_hat)
grad.x_hat=1i.*Kx.*func_hat;
grad.y_hat=1i.*Ky.*func_hat;