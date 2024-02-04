function [etainit,phiinit]=funG_IVP_data2d(data,dom)

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
else
phiinit=zeros(dom.Ny,dom.Nx);
end

end