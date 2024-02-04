function [Gen] = GenerationMethod0(g,Infmodel,dom,Om,Ug,ugadd,depth,...
    Spatinflux,Influx,ramp2d,I,timeSimul,tapered)

 timesig=Influx.wave(I).time;
 ww =funC_freqspace(timesig);

 XYinfline=Spatinflux.line(I).xy;
 ymin=dom.Y(funC_closest(dom.Y,min(XYinfline(:,2))));
 ymax=dom.Y(funC_closest(dom.Y,max(XYinfline(:,2))));
 xmin=dom.X(funC_closest(dom.X,min(XYinfline(:,1))));
 xmax=dom.X(funC_closest(dom.X,max(XYinfline(:,1))));
 
if strcmp(Spatinflux.line(I).Orientation,'Vertical')
    Insig=Influx.wave(I).eta;
    kk           = dom.kx;
    k1           = dom.ky;        

    Heav=funC_Heav2d(dom.Y,dom.fbl.b,dom.fbl.t,dom.Nx);
    Yinfb=ymin-dom.Y(1);
    Yinft=dom.Y(end)-ymax;
    Heav=Heav.*funC_Heav2d(dom.Y,Yinfb,Yinft,dom.Nx);
    if tapered.check==1
    lengthInf=ymax-ymin;
    taperlength.l=0;taperlength.r=0;
    taperlength.t=(round(tapered.length*lengthInf/dom.dy)+1)*dom.dy;
    taperlength.b=taperlength.t;

    taperChar=funC_tapered2d(dom.X,dom.Y,dom.X(1),dom.X(end),...
                  ymin,ymax,taperlength);               
    Heav=Heav.*taperChar;
    end
    
    

else
    Insig=Influx.wave(I).eta;
    kk            = dom.ky;
    k1            = dom.kx;
    Heav=funC_Heav2d(dom.X,dom.fbl.l,dom.fbl.r,dom.Ny);
    Xinfl=xmin-dom.X(1);Xinfr=dom.X(end)-xmax;
    Heav=Heav.*funC_Heav2d(dom.X,Xinfl,Xinfr,dom.Ny);
     if tapered.check==1
    lengthInf=xmax-xmin;
    taperlength.l=(round(tapered.length*lengthInf/dom.dx)+1)*dom.dx;
    taperlength.r=taperlength.l;
    taperlength.t=0;taperlength.b=0;
    taperChar=funC_tapered2d(dom.X,dom.Y,xmin,xmax,...
                  dom.Y(1),dom.Y(end),taperlength);  
    Heav=Heav.*taperChar';
    
     end
end

 K1=funOprt_invOmExact(k1,ww,depth,g);
[kk1,KK1]=meshgrid(k1,K1);
KK=sqrt(KK1.^2+kk1.^2);

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
        dxGam(indyinf,:)=funC_ifft(1i.*Om(kk,depth,g,ugadd).*fft(deltadir)); %for uni-direction
        end
        Gam_hat=fft2(Gam.*Heav)*fact;
        Gam      = funC_ifft2(Gam_hat);
        dxGam_hat=fft2(dxGam.*Heav)*fact;
        SkewGam    = funC_ifft2(dxGam_hat);
        
        GamKK1_hat=ones(size(KK1))/dx;
        %GamK1_hat=ones(size(K1))/dx;
    else
        Gam=zeros(dom.Ny,dom.Nx);dyGam=Gam;
        for jj=1:length(XYinfline(:,1))
        indXinf=funC_closest(dom.X,XYinfline(jj,1));
        indYinf=funC_closest(dom.Y,XYinfline(jj,2));
        deltadir=zeros(dom.Ny,1);
        deltadir(indYinf)=1/dy;
        Gam(:,indXinf)=deltadir;
      %  dyGam(:,indXinf)=funC_ifft(1i.*kk'.*fft(deltadir)); %for uni-direction
        dyGam(:,indXinf)=funC_ifft(1i.*Om(kk,depth,g,ugadd).*fft(deltadir)); %for uni-direction
       
        end
        Gam_hat =fft2(Gam.*Heav')*fact;
        Gam     = funC_ifft2(Gam_hat);
        dyGam_hat=fft2(dyGam.*Heav')*fact;    
        SkewGam    = funC_ifft2(dyGam_hat);
        
        GamKK1_hat=ones(size(KK1))/dy;
        %GamK1_hat=ones(size(K1))/dy;
    end
    
elseif strcmp(Infmodel.type(I),'Area')
    if strcmp(Spatinflux.line(I).Orientation,'Vertical')
        Gam=zeros(dom.Ny,dom.Nx);dxGam=Gam;
        for jj=1:length(XYinfline(:,1))
        indyinf=funC_closest(dom.Y,XYinfline(jj,2));
        Ugx=funC_ifft(Ug(kk,depth,g,ugadd));
        dist=XYinfline(jj,1)-dom.X(1);
        Ugx= circshift(Ugx,round(dist/(dx)));
        Gam(indyinf,:)      = Ugx;
        %dxGam(indyinf,:)=funC_ifft(1i.*kk'.*fft(Ugx)); %for uni-direction
        dxGam(indyinf,:)=funC_ifft(1i.*Om(kk,depth,g,ugadd).*fft(Ugx)); %for uni-direction
        end
        Gam_hat    = fft2(Gam.*Heav)*fact;
        Gam        = funC_ifft2(Gam_hat);
        %dxGam_hat  =1i.*funOprt_Om2dExact(dom.Kx,dom.Ky,depth,g,[0 1]).*Gam_hat;
        
        dxGam_hat  =fft2(dxGam.*Heav)*fact;
        SkewGam    = funC_ifft2(dxGam_hat);
        UgKK1      =funOprt_Ug2d(KK1,kk1,depth,g,[0 1]);
        UgxKK1     =funC_ifft2(UgKK1);
        GamKK1_hat = fft2(UgxKK1);
        
%         UgxK1=funC_ifft(Ug(K1,depth,ugadd));
%         GamK1_hat= fft(UgxK1);
    else
         Gam=zeros(dom.Ny,dom.Nx);dyGam=Gam;
        for jj=1:length(XYinfline(:,2))
        indxinf=funC_closest(dom.X,XYinfline(jj,1));
        Ugy      =funC_ifft(Ug(kk,depth,g,ugadd));
        dist     =XYinfline(jj,2)-dom.Y(1);
        Ugy      = circshift(Ugy,round(dist/(dy)));
        Gam(:,indxinf)= Ugy;
       % dyGam(:,indxinf)=funC_ifft(1i.*kk'.*fft(Ugy)); %for uni-direction
        dyGam(:,indxinf)=funC_ifft(1i.*Om(kk,depth,g,ugadd).*fft(Ugy)); %for uni-direction
        end
        Gam_hat  = fft2(Gam.*Heav')*fact;
        Gam      = funC_ifft2(Gam_hat);
        dyGam_hat= fft2(dyGam.*Heav')*fact;
        SkewGam    = funC_ifft2(dyGam_hat);
        
        UgKK1  =funOprt_Ug2d(kk1,KK1,depth,g,[1 0]);
        UgyKK1=funC_ifft2(UgKK1);
        GamKK1_hat= fft2(UgyKK1);
        
%         UgyK1=funC_ifft(Ug(K1,depth,ugadd));
%         GamK1_hat= fft(UgyK1);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(Spatinflux.line(I).Orientation,'Vertical')
UgKK_hat  =funOprt_Ug2d(KK1,kk1,depth,g,[0 1]);
ffact   =dom.dy./(2*pi);
else
UgKK_hat=funOprt_Ug2d(kk1,KK1,depth,g,[1 0]);  
ffact   =dom.dx./(2*pi);
end


%WW1=funOprt_Om2dExact(KK1,kk1,depth,g,[1 1]);
%Op_C=WW1./KK1;

Insig_hat=fft2(Insig)*ffact;
ff_hat=Insig_hat/(2*pi).*abs(KK1)./KK.*UgKK_hat./GamKK1_hat;
%Skewf_hat=Op_C.*ff_hat./(1i.*WW1);
ff=funC_ifft2(ff_hat);


Skewf=cumtrapz(timesig.',ff);
%Skewf=funC_ifft2(Skewf_hat);
f_time_spatial=ff.*ramp2d;
Skewf_time_spatial=Skewf.*ramp2d;
%%%%%%%%%%%%adjust to simulation time%%%%%%%%%%%%%%%%%%%%%%%%%
% if timeSimul.t_end>timesig(end)
% Nt0=length(timesig);
% dt=timesig(2)-timesig(1);
% Nadd=round((timeSimul.t_end-timesig(end))/dt)+1;
% timeSigtemp=zeros(Nt0+Nadd,1);
% f_time_spatialtemp=zeros(Nt0+Nadd,length(ff(1,:)));
% f_time_spatialtemp(1:Nt0,:)=f_time_spatial;
% f_time_spatial=f_time_spatialtemp;
% Skewf_time_spatial_temp=zeros(Nt0+Nadd,length(ff(1,:)));
% Skewf_time_spatial_temp(1:Nt0,:)=Skewf_time_spatial;
% Skewf_time_spatial=Skewf_time_spatial_temp;
% timeSigtemp(1:Nt0)=timesig;
% timeSigtemp(Nt0+1:end)=timesig(end)+dt.*(1:Nadd);
% timesig=timeSigtemp';
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Gen.Signal=spline(timesig,f_time_spatial');
Gen.Gam    =Gam;
Gen.SkewGam=SkewGam;
Gen.SkewSignal=spline(timesig,Skewf_time_spatial');


% if strcmp(Spatinflux.line(I).Orientation,'Vertical')
% [T,Y]=meshgrid(dom.Y,timesig);
% figure(180)
% subplot(2,1,1)
% surf(T,Y,Skewf,'edgecolor','none');%f_time_spatial'
% view(2);
% colorbar
% subplot(2,1,2)
% surf(T,Y,Skewf1,'edgecolor','none')
% view(2);
% colorbar
% else
% [X,T]=meshgrid(dom.X,timesig);
% figure(181)
% subplot(2,1,1)
% surf(X,T,Insig,'edgecolor','none')
% view(2);
% colorbar
% subplot(2,1,2)
% surf(X,T,f_time_spatial,'edgecolor','none')
% view(2);
% colorbar
% end

end







