function bathyshore = funB_Bathy_shore(dom,bath)
parslope = bath.par;

Depth       = parslope(1);               % deepest depth
sf          = parslope(2);              % steepness fraction


if sign(sf)==1
    Xse         = parslope(3);              % shore
    length_s    = abs(Depth)/abs(sf);            % length of slope
    Xss         = Xse - length_s;  % start slope
    
    D_left=Depth;D_right=0;flagL=1;flagR=0;
    if strcmp(bath.type,'Shore (in x-axis)')%%if strcmp(bath.type,'Sx') %%nida
        fbl.l = 0; fbl.r   = dom.fbl.r;
        fbl.t = 0; fbl.b   = 0;
        cfSA  = funC_cfSA2d(dom.X,dom.Y,fbl);
        Heav  =funC_Heaviside(dom.X-Xss);
    else
        fbl.l = 0;         fbl.r   = 0;
        fbl.t = dom.fbl.t; fbl.b   = 0;
        cfSA  = funC_cfSA2d(dom.X,dom.Y,fbl);
        Heav  =funC_Heaviside(dom.Y-Xss);
    end
else
    Xss         = parslope(3);              % shore
    length_s    = abs(Depth)/abs(sf);            % length of slope
    Xse         = Xss + length_s;  % start slope
    
    D_left=0;D_right=Depth;flagL=0;flagR=1;
    if strcmp(bath.type,'Shore (in x-axis)') %%if strcmp(bath.type,'Sx') %%nida
        fbl.l = dom.fbl.l; fbl.r   = 0;
        fbl.t = 0; fbl.b   = 0;
        cfSA  = funC_cfSA2d(dom.X,dom.Y,fbl);
        Heav  =(1-funC_Heaviside(dom.X-Xse));
%         figure;
%         plot(dom.X,cfSA,'r',dom.X,Heav,'b')
    else
        fbl.l = 0;  fbl.r   = 0;
        fbl.t = 0;  fbl.b   = dom.fbl.b;
        cfSA  = funC_cfSA2d(dom.X,dom.Y,fbl);
        Heav  =1-funC_Heaviside(dom.Y-Xse);
%         figure;
%         plot(dom.Y,cfSA,'r',dom.Y,Heav,'b')
    end
end

sf=abs(sf);
if strcmp(bath.type,'Shore (in x-axis)')%%if strcmp(bath.type,'Shore (in x-axes)')    %%nida
    bathyshore = flagL*D_left*ones(1,dom.Nx) + ...
        -sign(D_left-D_right)*sf*(dom.X-Xss).*Heav+ ...
        + flagR*D_right.*funC_Heaviside(dom.X-Xse);
           
%     figure;
%         plot(dom.X,bathyshore,'r')

    bathyshore =-repmat(bathyshore,dom.Ny,1).*cfSA;
%     figure;
%         surf(bathyshore,'edgecolor','none')
    
else
    bathyshore = flagL*D_left*ones(1,dom.Ny) + ...
        -sign(D_left-D_right)*sf*(dom.Y-Xss).*Heav + ...
        + flagR*D_right.*funC_Heaviside(dom.Y-Xse);
    bathyshore =-repmat(bathyshore',1,dom.Nx).*cfSA;
end
% 
% [X,Y]=meshgrid(dom.X,dom.Y);
% figure
% surf(X,Y,bathyshore,'edgecolor','none');