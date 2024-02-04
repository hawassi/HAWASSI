function [Gen] = GenerationMethod_bathy(g,Infmodel,dom,Ug,ugadd,depthInfl,...
    Spatinflux,Influx,ramp2d,I,InterpVg,timeSimul,tapered)

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
    
    Depth_influx_line=zeros(dom.Ny,1);
    Indymin=funC_closest(dom.Y,ymin);Indymax=funC_closest(dom.Y,ymax);
    Depth_influx_line(Indymin:Indymax)=depthInfl;
    Depth_influx_line(1:Indymin-1)=depthInfl(1);
    Depth_influx_line(Indymax+1:end)=depthInfl(end);
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
    Depth_influx_line=zeros(dom.Nx,1);
    Indxmin=funC_closest(dom.X,xmin);Indxmax=funC_closest(dom.X,xmax);
    Depth_influx_line(Indxmin:Indxmax)=depthInfl;
    Depth_influx_line(1:Indxmin-1)=depthInfl(1);
    Depth_influx_line(Indxmax+1:end)=depthInfl(end);
end

KK1=funOprt_invOmExact(k1,ww,Depth_influx_line,g);
kk1=repmat(k1,length(ww),1);
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
        dxGam(indyinf,:)=funC_ifft(1i.*kk'.*fft(deltadir)); %for uni-direction
        end
        Gam_hat=fft2(Gam.*Heav)*fact;
        Gam      = funC_ifft2(Gam_hat);
        dxGam_hat=fft2(dxGam'.*Heav)*fact;
        SkewGam    = funC_ifft2(dxGam_hat);
        
        GamKK1_hat=ones(size(KK1))/dx;
        %GamK1_hat=ones(size(K1))/dx;
    else
       Gam=zeros(dom.Ny,dom.Nx);dyGam=Gam;
        for jj=1:length(XYinfline(:,1))
        indXinf=funC_closest(dom.X,XYinfline(:,1));
        indYinf=funC_closest(dom.Y,XYinfline(:,2));
        deltadir=zeros(dom.Ny,1);
        deltadir(indYinf)=1/dy;
        Gam(:,indXinf)=deltadir;
        dyGam(:,indXinf)=funC_ifft(1i.*kk'.*fft(deltadir)); %for uni-direction
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
        Ugx=funC_ifft(Ug(kk,depthInfl(jj),g,ugadd));
        dist=XYinfline(jj,1)-dom.X(1);
        Ugx= circshift(Ugx,round(dist/(dx)));
        Gam(indyinf,:)      = Ugx;
        dxGam(indyinf,:)=funC_ifft(1i.*kk'.*fft(Ugx)); %for uni-direction
        end
        Gam_hat    = fft2(Gam.*Heav)*fact;
        Gam        = funC_ifft2(Gam_hat);
        dxGam_hat  =fft2(dxGam.*Heav)*fact;
        SkewGam    = funC_ifft2(dxGam_hat);
        
        %         UgKK1    =Oprt_Ug2d(KK1,kk1,depth,[0 1]);
        %         UgxKK1=funC_ifft2(UgKK1);
        %         GamKK1_hat= fft2(UgxKK1);
        
        %         UgxK1=funC_ifft(Ug(K1,depth,ugadd));
        %         GamK1_hat= fft(UgxK1);
    else
        Gam=zeros(dom.Ny,dom.Nx);dyGam=Gam;
        for jj=1:length(XYinfline(:,2))
        indxinf=funC_closest(dom.X,XYinfline(jj,1));
        Ugy      =funC_ifft(Ug(kk,depthInfl(jj),g,ugadd));
        dist     =XYinfline(jj,2)-dom.Y(1);
        Ugy      = circshift(Ugy,round(dist/(dy)));
        Gam(:,indxinf)= Ugy;
        dyGam(:,indxinf)=funC_ifft(1i.*kk'.*fft(Ugy)); %for uni-direction
        end
        Gam_hat  = fft2(Gam.*Heav')*fact;
        Gam      = funC_ifft2(Gam_hat);
        dyGam_hat= fft2(dyGam.*Heav')*fact;
        SkewGam    = funC_ifft2(dyGam_hat);
        
        %         UgKK1  =Oprt_Ug2d(kk1,KK1,depth,[1 0]);
        %         UgyKK1=funC_ifft2(UgKK1);
        %         GamKK1_hat= fft2(UgyKK1);
        
        %         UgyK1=funC_ifft(Ug(K1,depth,ugadd));
        %         GamK1_hat= fft(UgyK1);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(Infmodel.type(I),'Point')
    %%% Parameters for interpolating an operator
    if strcmp(Spatinflux.line(I).Orientation,'Vertical')
        Depth_influx_line=repmat(Depth_influx_line',length(ww),1);
    else
        Depth_influx_line=repmat(Depth_influx_line,length(ww),1);
    end
    gam_plus=interp1(InterpVg.DDepth,InterpVg.Ggam_plus,Depth_influx_line,'spline');
    gam_min=interp1(InterpVg.DDepth,InterpVg.Ggam_min,Depth_influx_line,'spline');
    gam_mid=interp1(InterpVg.DDepth,InterpVg.Ggam_mid,Depth_influx_line,'spline');
    Dmin=InterpVg.D_min;Dplus=InterpVg.D_plus;
    Dmid=InterpVg.D_mid;
    %%%
end
if strcmp(Spatinflux.line(I).Orientation,'Vertical')
    if strcmp(Infmodel.type(I),'Point')
        UgKK_hat_min=Oprt_Ug2d(KK1,kk1,Dmin,g,[0 1]);
        UgKK_hat_mid=Oprt_Ug2d(KK1,kk1,Dmid,g,[0 1]);
        UgKK_hat_plus=Oprt_Ug2d(KK1,kk1,Dplus,g,[0 1]);
    end
    ffact   =dom.dy./(2*pi);
else
    if strcmp(Infmodel.type(I),'Point')
        UgKK_hat_min=Oprt_Ug2d(kk1,KK1,Dmin,g,[1 0]);
        UgKK_hat_mid=Oprt_Ug2d(kk1,KK1,Dmid,g,[1 0]);
        UgKK_hat_plus=Oprt_Ug2d(kk1,KK1,Dplus,g,[1 0]);
    end
    ffact   =dom.dx./(2*pi);
end

if strcmp(Infmodel.type(I),'Point')
    Insig_hat=fft2(Insig)*ffact;
    ff_hat=Insig_hat/(2*pi).*abs(KK1)./KK./GamKK1_hat;
    ff    =funC_ifft2(ff_hat);
    ff_hat=funOprt_HSS2d(UgKK_hat_min,gam_min,ff_hat,ff)+...
        funOprt_HSS2d(UgKK_hat_plus,gam_plus,ff_hat,ff)+...
        funOprt_HSS2d(UgKK_hat_mid,gam_mid,ff_hat,ff);
else
    Insig_hat=fft2(Insig)*ffact;
    ff_hat=Insig_hat/(2*pi).*abs(KK1)./KK;
end
Skewf_hat=ff_hat./(1i.*KK1);

ff=funC_ifft2(ff_hat);
Skewf=funC_ifft2(Skewf_hat);
f_time_spatial=ff.*ramp2d;
Skewf_time_spatial=Skewf.*ramp2d;



%%%%%%%%%%%%additional time simulation%%%%%%%%%%%%%%%%%%%%%%%%%
% indt1=funC_closest(timesig,timeSimul.t_init);
% indt2=funC_closest(timesig,timeSimul.t_end);
% dtSig=timesig(2)-timesig(1);
% if timeSimul.dt>dtSig
%     tcoarse     = round(timeSimul.dt/dtSig);
% else tcoarse=1;
% end
% timesig=timesig(indt1:tcoarse:indt2);
% 
% if timeSimul.t_end>timesig(end)
%     Nt0=length(timesig);
%     dt=timesig(2)-timesig(1);
%     Nadd=round((timeSimul.t_end-timesig(end))/dt)+1;
%     timeSigtemp=zeros(Nt0+Nadd,1);
%     f_time_spatialtemp=zeros(Nt0+Nadd,length(ff(1,:)));
%     f_time_spatialtemp(1:Nt0,:)=f_time_spatial;
%     f_time_spatial=f_time_spatialtemp;
%     Skewf_time_spatial_temp=zeros(Nt0+Nadd,length(ff(1,:)));
%     Skewf_time_spatial_temp(1:Nt0,:)=Skewf_time_spatial;
%     Skewf_time_spatial=Skewf_time_spatial_temp;
%     timeSigtemp(1:Nt0)=timesig;
%     timeSigtemp(Nt0+1:end)=timesig(end)+dt.*(1:Nadd);
%     timesig=timeSigtemp';
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gen.Signal=spline(timesig,f_time_spatial');
Gen.Gam    =Gam;
Gen.SkewGam=SkewGam;
Gen.SkewSignal=spline(timesig,Skewf_time_spatial');

% if strcmp(Spatinflux.line(I).Orientation,'Vertical')
%     [T,Y]=meshgrid(timesig,dom.Y);
%     figure
%     subplot(2,1,1)
%     surf(T,Y,f_time_spatial','edgecolor','none')
%     view(2);
%     colorbar
%     subplot(2,1,2)
%     surf(T,Y,Skewf_time_spatial','edgecolor','none')
%     view(2);
%     colorbar
% else
%     [X,T]=meshgrid(dom.X,timesig);
%     figure
%     subplot(2,1,1)
%     surf(X,T,f_time_spatial,'edgecolor','none')
%     view(2);
%     colorbar
%     subplot(2,1,2)
%     surf(X,T,Skewf_time_spatial,'edgecolor','none')
%     view(2);
%     colorbar
% end

end







