function phitan=TangentialPotential2dbyMinimization1(g,dom,Oprt,chiB,phiin,phiinhat,zeta)

nu         = 0;
phitan     =phiin.*(1-chiB);
%phitan=TangentialPotential2d(dom,Oprt,chiB,phiin,phiinhat);
Error_fun=@(nu,DtNphi,chiB) trapz(dom.Y,trapz(dom.X,(DtNphi-nu).*chiB,2));
grad_Obj_fun=@(nu,DtNphi,chiB) (DtNphi-nu).*chiB;
tol=0.001;
error=1;
dtf=0.01;

while error>tol
    DtNphitan= funC_ifft2(funOprt_L2d(g,Oprt.Om2dSqInterp,fft2(phitan),phitan));
    %DtNphitan= funC_ifft2(funOprt_L2dAB2F(dom,Oprt,zeta,fft2(phitan)));
    gradF=grad_Obj_fun(nu,DtNphitan,chiB);
    phitan=phitan-gradF*dtf;
    DtNphitan= funC_ifft2(funOprt_L2d(g,Oprt.Om2dSqInterp,fft2(phitan),phitan));
    %DtNphitan= funC_ifft2(funOprt_L2dAB2F(dom,Oprt,zeta,fft2(phitan)));
    error=abs(Error_fun(nu,DtNphitan,chiB));
end


% StepS=4;
% [XX,YY]=meshgrid(dom.X(1:StepS:end),dom.Y(1:StepS:end));
% subplot(3,1,2)
% surf(XX,YY,phitan(1:StepS:end,1:StepS:end),'edgecolor','none');
% view(2)
% subplot(3,1,3)
% surf(XX,YY,DtNphitan(1:StepS:end,1:StepS:end),'edgecolor','none');
% view(2)