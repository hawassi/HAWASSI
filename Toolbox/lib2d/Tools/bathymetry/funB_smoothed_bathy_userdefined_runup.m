function bathy=funB_smoothed_bathy_userdefined_runup(dom,shoretype)
bathyCharplus=zeros(size(dom.XX));
bathyCharplus(dom.bathy.profile>0)=ones;

if strcmpi(shoretype,'Userdefined (shore in x-axis)')
    fbl.l = dom.fbl.l; fbl.r   = dom.fbl.r;
    fbl.t = 0; fbl.b   = 0;
else
    fbl.t = dom.fbl.t; fbl.b   = dom.fbl.b;
    fbl.l = 0; fbl.r   = 0;
end
cfSA  = funC_cfSA2d(dom.X,dom.Y,fbl);
cfSA(bathyCharplus==0)=1;
bathy=dom.bathy.profile.*cfSA;
% figure;
% plot(dom.X,bathy(1,:))
end