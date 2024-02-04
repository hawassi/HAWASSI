function [radpot,AddedMass]=RadiationPotential2dbyMinimization(g,dom,Oprt,body)
chiB=body.par.chiB(end).char;

switch body.RBdyn
    case 'heave'
        nu = chiB;
    case 'surge'
       
    case 'pitch'
end

radpot     =zeros(size(chiB));%.*chiB;

Error_fun=@(nu,DtNpsi,chiB) trapz(dom.Y,trapz(dom.X,(DtNpsi-nu).*chiB,2));
grad_Obj_fun=@(nu,DtNpsi,chiB) (DtNpsi-nu).*chiB;
tol=0.001;
error=1;
dtf=0.01;



while error>tol
    DtNradpot= funC_ifft2(funOprt_L2d(g,Oprt.Om2dSqInterp,fft2(radpot),radpot));
    gradF=grad_Obj_fun(nu,DtNradpot,chiB);
    radpot=radpot-gradF*dtf;
    DtNradpot= funC_ifft2(funOprt_L2d(g,Oprt.Om2dSqInterp,fft2(radpot),radpot));
    error=abs(Error_fun(nu,DtNradpot,chiB));
end

AddedMass=zeros(body.N,1);
for i=1:body.N
chiBI=body.par.chiB(i).char;
AddedMass(i)=trapz(dom.Y,trapz(dom.X,radpot.*chiBI,2));
end

     StepS=2;
    [XX,YY]=meshgrid(dom.X(1:StepS:end),dom.Y(1:StepS:end));
    figure;
    subplot(2,1,1)
    surf(XX,YY,radpot(1:StepS:end,1:StepS:end),'edgecolor','none');
    title('Radiation Potential')
    subplot(2,1,2)
    surf(XX,YY,DtNradpot(1:StepS:end,1:StepS:end),'edgecolor','none');
    title('DtN Radiation Potential')
    pause(0.01)

