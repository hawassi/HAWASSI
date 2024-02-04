function phitan=TangentialPotential2dbyMinimization(dom,Oprt,chiB,phiin,phiinhat)

DNphiin_hat=funOprt_L2d(Oprt.Om2dSqInterp,phiinhat,phiin);
DNphiin    =funC_ifft2(DNphiin_hat);
nu         = -DNphiin.*chiB;
%phidif_hat =funOprt_OneperL2d(Oprt.OneperOm2dSqInterp,fft2(nu),nu);
phidif     =zeros(size(chiB));%funC_ifft2(phidif_hat);%.*chiB;
% phidif_hat =fft2(phidif);
Error_fun=@(nu,DtNphi,chiB) trapz(dom.Y,trapz(dom.X,(DtNphi-nu).*chiB,2));
grad_Obj_fun=@(nu,DtNphi,chiB) (DtNphi-nu).*chiB;
tol=0.0001;
error=1;
dtf=0.01;
while error>tol
    DtNphidif= funC_ifft2(funOprt_L2d(Oprt.Om2dSqInterp,fft2(phidif),phidif));
    gradF=grad_Obj_fun(nu,DtNphidif,chiB);
    phidif=phidif-gradF*dtf;
    DtNphidif= funC_ifft2(funOprt_L2d(Oprt.Om2dSqInterp,fft2(phidif),phidif));
    error=abs(Error_fun(nu,DtNphidif,chiB));
end
phitan     =phiin+phidif.*chiB;

% DtNphitan= funC_ifft2(funOprt_L2d(Oprt.Om2dSqInterp,fft2(phitan),phitan));
% StepS=4;
% [XX,YY]=meshgrid(dom.X(1:StepS:end),dom.Y(1:StepS:end));
% subplot(3,1,2)
% surf(XX,YY,phitan(1:StepS:end,1:StepS:end),'edgecolor','none');
% view(2)
% subplot(3,1,3)
% surf(XX,YY,DtNphitan(1:StepS:end,1:StepS:end),'edgecolor','none');
% view(2)
