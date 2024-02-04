function INsig=funG_harmonic_oblique(g,ampl,Tp,dom,time,depth,Spatinflux_type,theta0)


W0=2*pi./Tp;

if max(depth)-min(depth)<1e-3
Depth=max(depth);
else
Depth=depth;
end


if strcmp(Spatinflux_type,'Vertical')
    K0=funOprt_invOmExact(dom.ky,W0,Depth,g);
    if length(K0)==dom.Ny
        K0=repmat(K0,length(time),1);
    end
   
    [Y,T]=meshgrid(dom.Y,time);
    INsig=ampl.*cos(K0.*sin(theta0).*Y-W0.*T);
    
else
    K0=funOprt_invOmExact(dom.kx,W0,Depth,g);
    if length(K0)==dom.Nx
        K0=repmat(K0,length(time),1);
    end
    [X,T]=meshgrid(dom.X,time);
    INsig=ampl.*cos(K0.*cos(theta0).*X-W0.*T);
end

% [X,Y]=meshgrid(dom.X,dom.Y);
% etaa=ampl.*cos(K0.*cos(theta0).*X+K0.*sin(theta0).*Y-W0.*20);

% figure;
% surf(X,Y,etaa,'edgecolor','none')
% xlabel('X');
% ylabel('Y');
% view(2);
% colorbar;