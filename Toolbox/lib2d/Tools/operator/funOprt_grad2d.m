function grad=funOprt_grad2d(Kx,Ky,func_hat)
grad.x=funC_ifft2(1i.*Kx.*func_hat);
grad.y=funC_ifft2(1i.*Ky.*func_hat);