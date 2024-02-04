function INsig=funG_oblique(ampl,Tp,dom,ww,depth,Spatinflux_type,theta0)


dww=ww(2)-ww(1);Nw=length(ww);

K1=invOmExact(ww,depth);

if strcmp(Spatinflux_type,'Vertical')
    dky=dom.ky(2)-dom.ky(1);
    [kk1,KK1]=meshgrid(dom.ky,K1);
    KK=sqrt(KK1.^2+kk1.^2);
    
    deltadir=zeros(size(KK));
    
    for ii=1:Nw
        K1wtheta=K1(ii).*sin(theta0);
        indky0=funC_closest(dom.ky,K1wtheta);
        deltadir(ii,indky0)=1/dky;
    end
    Nspat=dom.Ny;dk=dky;
else
    dkx=dom.kx(2)-dom.kx(1);
    [kk1,KK1]=meshgrid(dom.kx,K1);
    KK=sqrt(KK1.^2+kk1.^2);
    
    deltadir=zeros(size(KK));
    for ii=1:Nw
        K1wtheta=K1(ii).*sin(theta0);
        indkx0=funC_closest(dom.kx,K1wtheta);
        deltadir(ii,indkx0)=1/dkx;
    end
    Nspat=dom.Nx;dk=dkx;
end
W0=2*pi./Tp;
delta_W0=zeros(size(ww));
indW0=funC_closest(ww,W0);
delta_W0(indW0)=1./dww;
INsig_hat=zeros(size(deltadir));
for jj=1:dom.Ny
    INsig_hat(:,jj)=ampl.*deltadir(:,jj).*delta_W0;
end
INsig=funC_ifft2(INsig_hat)*Nspat.*Nw*dk*dww/2;


figure;
surf(deltadir,'edgecolor','none')
view(2);
colorbar;