function [Gam, SkewGam] = funW_FluxSpatGam2D_1(g,dom,Ug,ugadd)

bdyInflidx=dom.wall.bdyInfl.index;
bdyInfldirflag=dom.wall.bdyInfl.dirflag;
Gam=zeros(dom.Ny,dom.Nx);
SkewGam_min2=Gam;SkewGam_min1=Gam;
SkewGam_plus2=Gam;SkewGam_plus1=Gam;

dx  =dom.X(2)-dom.X(1);
dy  =dom.Y(2)-dom.Y(1);
fact=(2*pi).^2./(dx.*dy);


for ii=1:length(bdyInflidx)
    depth=-dom.bathy.profile(bdyInflidx(ii));
    if abs(bdyInfldirflag(ii))==1  %%% in x axis
        kk           = dom.kx;
        Ugx=funC_ifft(Ug(kk,depth,g,ugadd));
        indyinf=funC_closest(dom.Y,dom.YY(bdyInflidx(ii)));
        dist   = dom.XX(bdyInflidx(ii))-dom.X(1);
        Ugx    = circshift(Ugx,round(dist/(dx)));
        Gam(indyinf,:)    = Gam(indyinf,:)+Ugx'.*fact;
        tempSkew          =-sign(bdyInfldirflag(ii)).*funC_ifft(1i.*OmExact(kk,depth).*fft(Ugx)).'.*fact;
        if bdyInfldirflag(ii)==-1
            SkewGam_min1(indyinf,:)= SkewGam_min1(indyinf,:)+tempSkew;
        else
            SkewGam_plus1(indyinf,:)= SkewGam_plus1(indyinf,:)+tempSkew;
        end
    elseif abs(bdyInfldirflag(ii))==2
        kk            = dom.ky; %%% in y axis
        indxinf=funC_closest(dom.X,dom.XX(bdyInflidx(ii)));
        Ugy      =funC_ifft(Ug(kk,depth,g,ugadd));
        dist     =dom.YY(bdyInflidx(ii))-dom.Y(1);
        Ugy      = circshift(Ugy,round(dist/(dy)));
        Gam(:,indxinf)=Gam(:,indxinf)+ Ugy.*fact;
        tempSkew     =-sign(bdyInfldirflag(ii)).*funC_ifft(1i.*OmExact(kk,depth).*fft(Ugy)).*fact; %for uni-direction
        if bdyInfldirflag(ii)==-2
            SkewGam_min2(:,indxinf)= SkewGam_min2(:,indxinf)+tempSkew;
        else
            SkewGam_plus2(:,indxinf)= SkewGam_plus2(:,indxinf)+tempSkew;
        end
    end
end
Gam=Gam.*(dom.wall.char);

SkewGam.min2=SkewGam_min2.*(dom.wall.char);
SkewGam.min1=SkewGam_min1.*(dom.wall.char);
SkewGam.plus1=SkewGam_plus1.*(dom.wall.char);
SkewGam.plus2=SkewGam_plus2.*(dom.wall.char);
SkewGam.All  =SkewGam.min2+SkewGam.min1+SkewGam.plus1+SkewGam.plus2;
end







