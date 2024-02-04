function [etainit,uinit,vinit]=funG_IVP_data2duv(data,dom)
uinit=zeros(dom.Ny,dom.Nx);vinit=zeros(dom.Ny,dom.Nx);

Xdat=data(:,1);Ydat=data(:,2);Etadat=data(:,3);

xx=reshape(dom.XX,1,[]);
yy=reshape(dom.YY,1,[]);
Feta=scatteredInterpolant(Xdat,Ydat,Etadat,'nearest');
Eta=Feta(xx,yy);
etainit=reshape(Eta,dom.Ny,[]).*dom.cfSA;
if length(data(1,:))==4
    Phidat=data(:,4);
    Fphi=scatteredInterpolant(Xdat,Ydat,Phidat,'nearest');
    Phi=Fphi(xx,yy);
    phiinit=reshape(Phi,dom.Ny,[]).*dom.cfSA;
    uinit=funC_ifft2(1i.*dom.Kx.*fft2(phiinit)).*dom.cfSA;
    vinit=funC_ifft2(1i.*dom.Ky.*fft2(phiinit)).*dom.cfSA;
elseif length(data(1,:))==5
    udat=data(:,4);
    Fu=scatteredInterpolant(Xdat,Ydat,udat,'nearest');
    u=Fu(xx,yy);
    uinit=reshape(u,dom.Ny,[]).*dom.cfSA;
    vdat=data(:,5);
    Fv=scatteredInterpolant(Xdat,Ydat,vdat,'nearest');
    v=Fv(xx,yy);
    vinit=reshape(v,dom.Ny,[]).*dom.cfSA;
end

end