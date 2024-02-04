function phitan=TangentialPotential2d(dom,Oprt,chiB,phiin,phiinhat)

DNphiin_hat= funOprt_L2d(Oprt.Om2dSqInterp,phiinhat,phiin);
DNphiin    = funC_ifft2(DNphiin_hat);
nu         = -DNphiin.*chiB;
phidif_hat = funOprt_OneperL2d(Oprt.OneperOm2dSqInterp,fft2(nu),nu);
phitan     = phiin+funC_ifft2(phidif_hat).*chiB;

% StepS=4;
% [XX,YY]=meshgrid(dom.X(1:StepS:end),dom.Y(1:StepS:end));
% subplot(3,1,2)
% surf(XX,YY,phiin(1:StepS:end,1:StepS:end),'edgecolor','none');
% subplot(3,1,3)
% surf(XX,YY,phitan(1:StepS:end,1:StepS:end),'edgecolor','none');
