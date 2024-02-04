function al=funSP_alias2d(dom)
kx=dom.kx;ky=dom.ky;
IdkxIn=floor(dom.Nx/dom.cutfracx);
IdkyIn=floor(dom.Ny/dom.cutfracy);
kxmax=kx(max(IdkxIn,1));
kymax=ky(max(IdkyIn,1));

[kX,kY]=meshgrid(kx,ky);

if IdkxIn>1 && IdkyIn>1
    Kmax=sqrt(kxmax^2+kymax^2);
    al=ones(size(kX));
    K=sqrt(kX.^2+kY.^2);
    al(K>Kmax)=zeros;
else
    alx=ones(size(kX));aly=alx;
    alx(abs(kX)>kxmax)=0;
    aly(abs(kY)>kymax)=0;
    al=alx.*aly;   
end


% figure;
% %subplot(2,1,1)
% surf(kX,kY,al,'edgecolor','none');view(2);
% % subplot(2,1,2)
% % surf(kX,kY,alSq,'edgecolor','none');view(2);


% 
% 
% 
