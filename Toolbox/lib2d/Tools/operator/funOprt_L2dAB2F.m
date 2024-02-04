function L2dAB2_Z_hat=funOprt_L2dAB2F(g,dom,Oprt,zeta,Z_hat)
gradZ_hat          = funOprt_grad2d_hat(dom.Kx,dom.Ky,Z_hat);
Cp2gradzetaL0Z_hat.x_hat = Oprt.Cp2dSq.*gradZ_hat.x_hat;
Cp2gradzetaL0Z_hat.y_hat = Oprt.Cp2dSq.*gradZ_hat.y_hat;
L0Z_hat         = -funOprt_div2d(dom.Kx,dom.Ky,Cp2gradzetaL0Z_hat)/g;
L0Z             = funC_ifft2(L0Z_hat);
gradzetaL0Z_hat = funOprt_grad2d_hat(dom.Kx,dom.Ky,fft2(zeta.*L0Z));
Cp2gradzetaL0Z_hat.x_hat = Oprt.Cp2dSq.*gradzetaL0Z_hat.x_hat;
Cp2gradzetaL0Z_hat.y_hat = Oprt.Cp2dSq.*gradzetaL0Z_hat.y_hat;

L0zetaLZ_hat  = -funOprt_div2d(dom.Kx,dom.Ky,Cp2gradzetaL0Z_hat)/g;

gradZ         = funOprt_grad2d(dom.Kx,dom.Ky,Z_hat);
zetagradphi_hat.x_hat=fft2(zeta.*gradZ.x);
zetagradphi_hat.y_hat=fft2(zeta.*gradZ.y);
divzetagradphi_hat = funOprt_div2d(dom.Kx,dom.Ky,zetagradphi_hat);
L2dAB2_Z_hat  =L0Z_hat-(L0zetaLZ_hat+divzetagradphi_hat);
% % if isnan(zeta)
% %    pause 
% % end
% L0Z_hat       =Oprt.L2d.*Z_hat;
% zetaL0Z_hat   = fft2(zeta.*funC_ifft2(L0Z_hat));
% L0zetaLZ_hat  = Oprt.L2d.*zetaL0Z_hat;
% gradZ       = funOprt_grad2d(dom.Kx,dom.Ky,Z_hat);
% zetagradZ_hat.x_hat=fft2(zeta.*gradZ.x);
% zetagradZ_hat.y_hat=fft2(zeta.*gradZ.y);
% divzetagradphi_hat = funOprt_div2d(dom.Kx,dom.Ky,zetagradZ_hat);
% L2dAB2_Z_hat  =L0Z_hat-(L0zetaLZ_hat+divzetagradphi_hat);
% % 
% plotVar1=funC_ifft2(L2dAB2_Z_hat);
% plotVar2=funC_ifft2(divzetagradphi_hat);
% StepS=4;
% [XX,YY]=meshgrid(dom.X(1:StepS:end),dom.Y(1:StepS:end));
% subplot(3,1,2)
% surf(XX,YY,plotVar1(1:StepS:end,1:StepS:end),'edgecolor','none');
% view(2)
% subplot(3,1,3)
% surf(XX,YY,plotVar2(1:StepS:end,1:StepS:end),'edgecolor','none');
% view(2)
