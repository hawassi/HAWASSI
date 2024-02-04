function bathy = funB_Bathy_userdefined(dom,bath)
Xdat=bath.userdata(:,1);Ydat=bath.userdata(:,2);
Bathdat=bath.userdata(:,3);
xx=reshape(dom.XX,1,[]);
yy=reshape(dom.YY,1,[]);
F=scatteredInterpolant(Xdat,Ydat,Bathdat,'nearest');
bathy=F(xx,yy);
bathy=reshape(bathy,dom.Ny,[]);

% w = warning('query','last');  %turn off warning 
% warning('off',w.identifier);
% bathy=griddata(Xdat,Ydat,Bathdat,dom.XX,dom.YY);
% if max(max(isnan(bathy)))==1
% w = warning('query','last');  %turn off warning 
% warning('off',w.identifier);
% bathy(isnan(bathy))=griddata(Xdat,Ydat,Bathdat,dom.XX(isnan(bathy)),dom.YY(isnan(bathy)),'nearest');
% end


% 
% hf0=figure;
% set(0,'CurrentFigure',hf0);
% set(gcf,'Renderer','zbuffer');
% surf(dom.XX,dom.YY,bathy,'edgecolor','none');
% view(2)
