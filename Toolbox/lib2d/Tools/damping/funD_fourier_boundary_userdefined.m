function cfSA=funD_fourier_boundary_userdefined(dom,userdata,smoothfact)
data.x=userdata(:,1);    
data.y=userdata(:,2);
 
mappedindex=funC_mappedmatrixindexes(dom.X,dom.Y,data.x,data.y);
cfSA=1-inpolygon(dom.XX,dom.YY,dom.XX(mappedindex),dom.YY(mappedindex));
cfSAA=zeros(1,length(dom.Y),length(dom.X));
cfSAA(1,:,:)=cfSA;

cfSAAsmooth=smooth3(cfSAA,'box',smoothfact);
cfSA=squeeze(cfSAAsmooth(1,:,:));

% hf0=figure(22);
% set(0,'CurrentFigure',hf0);
% set(gcf,'Renderer','zbuffer'); %due to graphics driver
% subplot(2,1,1)
% plot(dom.XX(mappedindex),dom.YY(mappedindex))
% subplot(2,1,2)
% surf(dom.XX,dom.Y,cfSA,'edgecolor','none');

