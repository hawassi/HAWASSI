function Eta=funBAssim_interpEtatoSimDom(xx,yy,bdyeta,dom)
[xx,yy]=meshgrid(xx,yy);
xx=reshape(xx,[],1);
yy=reshape(yy,[],1);
bdyeta=reshape(bdyeta,[],1);
F=scatteredInterpolant(xx,yy,bdyeta,'linear');
XX=reshape(dom.XX,[],1);
YY=reshape(dom.YY,[],1);
Eta=F(XX,YY);
Eta=reshape(Eta,dom.Ny,[]);
end