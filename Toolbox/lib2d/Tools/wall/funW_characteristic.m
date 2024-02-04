%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%    HAWASSI 2D                                          %%%%%%%%%%
%%%%%%%%%    Hamiltonian Wave-Ship-Structure Interaction         %%%%%%%%%%
%%%%%%%%%    copyright (c): LabMath-Indonesia                    %%%%%%%%%%
%%%%%%%%%    version 12 August 2016                              %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Characteristic function for wall                            %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [wallchar,gradwallchar]=funW_characteristic(dom,nonlin)
wallparam=dom.wall.param;
wallchar=1;
Heav=ones(dom.Ny,dom.Nx);
[XX,YY]=meshgrid(dom.X,dom.Y);
for I=1:length(wallparam(:,1))
if strcmpi(wallparam(I,1),'Rectangle')
bb=cell2mat(wallparam(I,7))-dom.Y(1);
tt=dom.Y(end)-cell2mat(wallparam(I,8));
HeavY=funC_Heav2d(dom.Y,bb,tt,dom.Nx);
ll=cell2mat(wallparam(I,5))-dom.X(1);
rr=dom.X(end)-cell2mat(wallparam(I,6));
HeavX=funC_Heav2d(dom.X,ll,rr,dom.Ny)';
Heav=1-HeavX.*HeavY;
elseif strcmpi(wallparam(I,1),'Circle')
R=cell2mat(wallparam(I,4));
xc=cell2mat(wallparam(I,2));
yc=cell2mat(wallparam(I,3));
Heav((XX-xc).^2+(YY-yc).^2<=R.^2)=0;
elseif strcmpi(wallparam(I,1),'User-defined')
walluserdata=dom.wall.userdata(I).bdry;
mappedindex=funC_mappedmatrixindexes(dom.X,dom.Y,walluserdata(:,1),walluserdata(:,2));
Heav=1-inpolygon(dom.XX,dom.YY,dom.XX(mappedindex),dom.YY(mappedindex));
end
Cutoff=1-funW_ReflCoef_to_cutoff(cell2mat(wallparam(I,10)),nonlin);
Heav(Heav==0)=Cutoff;
wallchar=(wallchar.*Heav);
[gradwallcharF.x,gradwallcharF.y]=funOprt_FDGradient(wallchar,dom.dx,dom.dy,1);
[gradwallcharB.x,gradwallcharB.y]=funOprt_FDGradient(wallchar,dom.dx,dom.dy,-1);
gradwallcharF.x(gradwallcharF.x>0)=0;gradwallcharF.y(gradwallcharF.y>0)=0;
gradwallcharB.x(gradwallcharB.x<0)=0;gradwallcharB.y(gradwallcharB.y<0)=0;
gradwallchar.x=(gradwallcharF.x+gradwallcharB.x);
gradwallchar.y=(gradwallcharF.y+gradwallcharB.y);
% [gradgradwallchar.x,gradgradwallchar.y]=funOprt_FDGradient(gradwallchar.x,dom.dx,dom.dy,0);
% [gradgradwallcharB.x,gradgradwallcharB.y]=funOprt_FDGradient(gradwallchar.x,dom.dx,dom.dy,-1);
% gradgradwallcharF.x(gradwallcharF.x>0)=0;gradwallcharF.y(gradwallcharF.y>0)=0;

% gradwallchar.x(find(gradwallchar.x,dom.Ny,'first'))=2*gradwallchar.x(find(gradwallchar.x,dom.Ny,'first'));
% gradwallchar.x(find(gradwallchar.x,dom.Ny,'last'))=0;

% gradwallchar.x=circshift(gradwallchar.x,[0 -1]);

% IndxWall=funC_indx_bdy(wallchar);
% gradwallchar.x=zeros(size(dom.XX));
% gradwallchar.x(IndxWall)=1./dom.dx.*(1-Cutoff);
% gradwallchar.y(IndxWall)=1./dom.dy.*(1-Cutoff);
end

% [X,Y]=meshgrid(dom.X,dom.Y);
% figure;
% set(gcf,'Renderer','zbuffer'); %due to graphics driver
% subplot(2,1,1)
% surf(X,Y,gradgradwallchar.x,'Edgecolor','none');
% shading interp;
% subplot(2,1,2)
% surf(X,Y,gradgradwallchar.y,'Edgecolor','none');
% shading interp;