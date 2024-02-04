%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%    HAWASSI 2D                                          %%%%%%%%%%
%%%%%%%%%    Hamiltonian Wave-Ship-Structure Interaction         %%%%%%%%%%
%%%%%%%%%    copyright (c): LabMath-Indonesia                    %%%%%%%%%%
%%%%%%%%%    version 12 August 2016                              %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Characteristic function for wall                            %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [wallcharAll,wallchar,gradwallchar,wallcharInfl,NwallInfl,wallReflCoef]=funW_characteristic_for_combined_method(dom,nonlin)
wallparam=dom.wall.param;
wallchar=ones(size(dom.XX));wallReflCoefUnif=0;wallcharInfl=ones(size(dom.XX));

FlagFreqDep=0;NwallInfl=0;
[XX,YY]=meshgrid(dom.X,dom.Y);
wallReflCoef.AllOneFlagCutK='True';
for I=1:length(wallparam(:,1))  
Heav=ones(dom.Ny,dom.Nx);    
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


if strcmpi(wallparam(I,9),'Influxing')
    NwallInfl=NwallInfl+1;
    if strcmpi(wallparam(I,10),'Uniform')
        wallReflCoefUnif=wallReflCoefUnif+(1-Heav).*cell2mat(wallparam(I,11));
        wallRelfCoefFreqDep(I).char=0;
        wallRelfCoefFreqDep(I).flag=0;
    else
        FlagFreqDep=1;
        wallRelfCoefFreqDep(I).char=1-Heav;
        wallRelfCoefFreqDep(I).flag=1;
    end
    wallcharInfl=wallcharInfl.*Heav;
else
    wallRelfCoefFreqDep(I).char=0;
    wallRelfCoefFreqDep(I).flag=0;
    
    Cutoff=1-funW_ReflCoef_to_cutoff(cell2mat(wallparam(I,11)),nonlin);
    Heav(Heav==0)=Cutoff;
    wallchar=(wallchar.*Heav);
    %%% restriction to energy truncation model for partial reflecting wall. If
    %%% any refl coef <1 energy trucation model uses a non expanded one
    %%% (DtN=-(C*(Cdx\phi.*Wallchar))/g), for fully reflecting wall uses an
    %%% expanded model (DtN=Lphi.*Wallchar-Flux.*delta);
    if cell2mat(wallparam(I,11))<1
        wallReflCoef.AllOneFlagCutK='False';
    end
end

  
end

[gradwallcharF.x,gradwallcharF.y]=funOprt_FDGradient(wallchar,dom.dx,dom.dy,1);
[gradwallcharB.x,gradwallcharB.y]=funOprt_FDGradient(wallchar,dom.dx,dom.dy,-1);
gradwallcharF.x(gradwallcharF.x>0)=0;gradwallcharF.y(gradwallcharF.y>0)=0;
gradwallcharB.x(gradwallcharB.x<0)=0;gradwallcharB.y(gradwallcharB.y<0)=0;
gradwallchar.x=(gradwallcharF.x+gradwallcharB.x);
gradwallchar.y=(gradwallcharF.y+gradwallcharB.y);

wallReflCoefUnif(wallReflCoefUnif>1)=1;


wallReflCoef.uniform=wallReflCoefUnif;
wallReflCoef.FlagFreqDep=FlagFreqDep;
wallReflCoef.FreqDep=wallRelfCoefFreqDep;

wallcharAll=1-((1-wallchar)+(1-wallcharInfl));

% [X,Y]=meshgrid(dom.X,dom.Y);
% ff=figure(100);
% set(ff,'Renderer','zbuffer'); %due to graphics driver
% subplot(3,1,1)
% surf(X,Y,wallchar,'Edgecolor','none');
% view(2)
% subplot(3,1,2)
% surf(X,Y,gradwallchar.x,'Edgecolor','none');
% view(2)
% subplot(3,1,3)
% surf(X,Y,gradwallchar.y,'Edgecolor','none');
% view(2)