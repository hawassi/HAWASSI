%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%    HAWASSI 2D                                          %%%%%%%%%%
%%%%%%%%%    Hamiltonian Wave-Ship-Structure Interaction         %%%%%%%%%%
%%%%%%%%%    copyright (c): LabMath-Indonesia                    %%%%%%%%%%
%%%%%%%%%    version 12 August 2016                              %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Characteristic function for wall                            %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function wallchar=funW_characteristic_shiftedOneGrid(dom,nonlin)
wallparam=dom.wall.param;
wallchar=1;
Heav=ones(dom.Ny,dom.Nx);
dx=dom.X(2)-dom.X(1);dy=dom.Y(2)-dom.Y(1);
[XX,YY]=meshgrid(dom.X,dom.Y);
for I=1:length(wallparam(:,1))
if strcmpi(wallparam(I,1),'Rectangle')
bb=cell2mat(wallparam(I,7))-dom.Y(1);
tt=dom.Y(end)-cell2mat(wallparam(I,8));
HeavY=funC_Heav2d(dom.Y,bb,tt,dom.Nx);
ll=cell2mat(wallparam(I,5))-dom.X(1)-dx;
rr=dom.X(end)-cell2mat(wallparam(I,6))+dx;
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
end

% [X,Y]=meshgrid(dom.X,dom.Y);
% figure;
% set(gcf,'Renderer','zbuffer'); %due to graphics driver
% surf(X,Y,Heav,'Edgecolor','none');
% shading interp;