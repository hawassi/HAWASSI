function spec=funSP_spec_from_images(eta,X,Y,depth,g)

    
dx                           =X(2)-X(1);
Nx                           =length(X);
kx                           =funC_freqspace(X);
kx1                          =ifftshift(kx);
if kx1(1)>kx1(2)
    kx1=[kx1(2:end);kx1(1)];
end
dy                           =Y(2)-Y(1);
Ny                           =length(Y);
ky                           =funC_freqspace(Y);
ky1                          =ifftshift(ky);
if ky1(1)>ky1(2)
    ky1=[ky1(2:end);ky1(1)];
end
dk                           =kx(2)-kx(1);
dtheta                       =deg2rad(0.1);
kmax                         =sqrt(max(kx)^2+max(ky)^2);
spec.k                            =0:dk:kmax;
spec.theta                        =0:dtheta:2*pi;
dw                                =2*pi/15000;
wmax                              =pi;
spec.w                            =0:dw:wmax;
[K,~]                         =meshgrid(spec.k,spec.theta);


fac                       =sqrt(2)/sqrt(Nx)/sqrt(Ny)/sqrt(max(kx))/sqrt(max(ky))/sqrt(pi);%sqrt(2*dx*dy/Nx/Ny/pi/pi/pi);%
spec.speckxky             =abs(fac*fft2(eta)).^2;
spec.specktheta           =K.*funC_CarttoPolCordinate(kx1,ky1,ifftshift(spec.speckxky),spec.k,spec.theta);
spec.w                    =funOprt_OmExact(spec.k,depth,g,[]);
spec.specwtheta           =zeros(length(spec.theta),length(spec.w));
Cg                        =funOprt_UgExact(spec.k,depth,g,[]);
for j=1:length(spec.w)
    spec.specwtheta(:,j)=spec.specktheta(:,j)./Cg(j);
    spec.specw(j)=trapz(spec.theta,spec.specwtheta(:,j));
end
% [~,ind]=max(spec.specw);
[~,ind]=max(smooth(spec.specw));
spec.Tp=2*pi/spec.w(ind);
spec.var=trapz(spec.w,spec.specw);
spec.Hs=4*sqrt(spec.var);
