function [Gam, SkewGam] = funW_FluxSpatGam2D(g,Infmodel,dom,Ug,ugadd,depth,...
    Spatinflux,I,tapered)

%  timesig=Influx.wave(I).time;
%  ww =funC_freqspace(timesig);

 XYinfline=Spatinflux.line(I).xy;
%  ymin=dom.Y(funC_closest(dom.Y,min(XYinfline(:,2))));
%  ymax=dom.Y(funC_closest(dom.Y,max(XYinfline(:,2))));
%  xmin=dom.X(funC_closest(dom.X,min(XYinfline(:,1))));
%  xmax=dom.X(funC_closest(dom.X,max(XYinfline(:,1))));
       
  Heav=1;     
if strcmp(Spatinflux.line(I).Orientation,'Vertical')
    kk           = dom.kx;
else
    kk            = dom.ky;
end

% K1=funOprt_invOmExact(k1,ww,depth,g);
% [kk1,KK1]=meshgrid(k1,K1);
% KK=sqrt(KK1.^2+kk1.^2);

dx  =dom.X(2)-dom.X(1);
dy  =dom.Y(2)-dom.Y(1);
fact=(2*pi)^2/(dx*dy);

if strcmp(Infmodel.type(I),'Point')
    
    if strcmp(Spatinflux.line(I).Orientation,'Vertical')
        Gam=zeros(dom.Ny,dom.Nx);dxGam=Gam;
        for jj=1:length(XYinfline(:,2))
        indyinf=funC_closest(dom.Y,XYinfline(jj,2));
        indXinf=funC_closest(dom.X,XYinfline(jj,1));
        deltadir=zeros(dom.Nx,1);
        deltadir(indXinf)=1/dx;
        Gam(indyinf,:)  =deltadir;
        dxGam(indyinf,:)=funC_ifft(1i.*OmExact(kk,depth).*fft(deltadir)); %for uni-direction
        end
        Gam_hat=fft2(Gam.*Heav)*fact;
        Gam      = funC_ifft2(Gam_hat);
        dxGam_hat=fft2(dxGam.*Heav)*fact;
        SkewGam    = funC_ifft2(dxGam_hat);
        
        %GamKK1_hat=ones(size(KK1))/dx;
        %GamK1_hat=ones(size(K1))/dx;
    else
        Gam=zeros(dom.Ny,dom.Nx);dyGam=Gam;
        for jj=1:length(XYinfline(:,1))
        indXinf=funC_closest(dom.X,XYinfline(jj,1));
        indYinf=funC_closest(dom.Y,XYinfline(jj,2));
        deltadir=zeros(dom.Ny,1);
        deltadir(indYinf)=1/dy;
        Gam(:,indXinf)=deltadir;
        dyGam(:,indXinf)=funC_ifft(1i.*OmExact(kk,depth).*fft(deltadir)); %for uni-direction
        end
        Gam_hat =fft2(Gam.*Heav')*fact;
        Gam     = funC_ifft2(Gam_hat);
        dyGam_hat=fft2(dyGam.*Heav')*fact;    
        SkewGam    = funC_ifft2(dyGam_hat);
        
        %GamKK1_hat=ones(size(KK1))/dy;
        %GamK1_hat=ones(size(K1))/dy;
    end
    
elseif strcmp(Infmodel.type(I),'Area')
    if strcmp(Spatinflux.line(I).Orientation,'Vertical')
        Gam=zeros(dom.Ny,dom.Nx);dxGam=Gam;
        for jj=1:length(XYinfline(:,1))
        indyinf=funC_closest(dom.Y,XYinfline(jj,2));
        Ugx=funC_ifft(Ug(kk,depth,g,ugadd));
        dist=XYinfline(jj,1)-dom.X(1);
        Ugx= circshift(Ugx,floor(dist/(dx)));
        Gam(indyinf,:)      = Ugx;
        dxGam(indyinf,:)=funC_ifft(1i.*OmExact(kk,depth).*fft(Ugx)); %for uni-direction
        end
        Gam_hat    = fft2(Gam.*Heav)*fact;
        Gam        = funC_ifft2(Gam_hat);
        dxGam_hat  =fft2(dxGam.*Heav)*fact;
        SkewGam    = funC_ifft2(dxGam_hat);
%         UgKK1      =funOprt_Ug2d(KK1,kk1,depth,g,[0 1]);
%         UgxKK1     =funC_ifft2(UgKK1);
        %GamKK1_hat = fft2(UgxKK1);
        
%         UgxK1=funC_ifft(Ug(K1,depth,ugadd));
%         GamK1_hat= fft(UgxK1);
    else
        disp ('Here')
         Gam=zeros(dom.Ny,dom.Nx);dyGam=Gam;
        for jj=1:length(XYinfline(:,2))
        indxinf=funC_closest(dom.X,XYinfline(jj,1));
        Ugy      =funC_ifft(Ug(kk,depth,g,ugadd));
        dist     =XYinfline(jj,2)-dom.Y(1);
        Ugy      = circshift(Ugy,floor(dist/(dy)));
        Gam(:,indxinf)= Ugy;
        dyGam(:,indxinf)=funC_ifft(1i.*OmExact(kk,depth).*fft(Ugy)); %for uni-direction
        end
        Gam_hat  = fft2(Gam.*Heav')*fact;
        Gam      = funC_ifft2(Gam_hat);
        dyGam_hat= fft2(dyGam.*Heav')*fact;
        SkewGam    = funC_ifft2(dyGam_hat);
        
%         UgKK1  =funOprt_Ug2d(kk1,KK1,depth,g,[1 0]);
%         UgyKK1=funC_ifft2(UgKK1);
       % GamKK1_hat= fft2(UgyKK1);
        
%         UgyK1=funC_ifft(Ug(K1,depth,ugadd));
%         GamK1_hat= fft(UgyK1);
    end
end






% if strcmp(Spatinflux.type(I),'Vertical')
% [T,Y]=meshgrid(Influx.timesig,dom.Y);
% figure
% subplot(2,1,1)
% surf(T,Y,f_time_spatial','edgecolor','none')
% view(2);
% colorbar
% subplot(2,1,2)
% surf(T,Y,Skewf_time_spatial','edgecolor','none')
% view(2);
% colorbar
% else
% [X,T]=meshgrid(dom.X,Influx.timesig);
% figure
% subplot(2,1,1)
% surf(X,T,f_time_spatial,'edgecolor','none')
% view(2);
% colorbar
% subplot(2,1,2)
% surf(X,T,Skewf_time_spatial,'edgecolor','none')
% view(2);
% colorbar
% end

end







