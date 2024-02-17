% 26 Oct 2015 ship shapes
function [parS,eta] = shapereformsetup0(x,eta,shippar,SXc,SZc)
% draft=shippar.draft+shippar.Z;
% ell=shippar.length/2;
Nship=shippar.Nship;
shapeX0Z0=shippar.form.shapeX0Z0;

shapeXcZc=zeros(Nship+1,length(x));
dxshapeXcZc=zeros(Nship+1,length(x));

nu.x=zeros(Nship,length(x)); % normal
nu.z=zeros(Nship,length(x)); % normal
nu.The=zeros(Nship,length(x)); % normal

dx=x(2)-x(1);
Chi=zeros(Nship+1,length(x));
ChiWl=zeros(Nship,length(x));
xShip  =zeros(Nship,6);
SChi=0;SshapeXcZc=0;SdxshapeXcZc=0;
SChiWl=0;
nuz=0;nux=0;
 for i=1:Nship
    shapeXcZ0         = circshift(shapeX0Z0(i,:),[floor(SXc(i)/dx),0]);
    shapeXcZc_shifted = (shapeXcZ0+SZc(i));
    Chi(i,:)    = Heaviside(-sign(shapeXcZ0(i,:))-10^(-12));
    shapeXcZc(i,:)=shapeXcZc_shifted.*Chi(i,:);
    dxshapeXcZc(i,:)= gradient(shapeXcZc(i,:),dx);
    indxL       = find(Chi(i,:),1,'first')-1;
  
    xShip(i,1)  = x(indxL);
    indxR       = find(Chi(i,:),1,'last')+1;
    xShip(i,3)  = x(indxR);
    indxC       = indxL+floor((indxR-indxL)/2);
    xShip(i,2)  = x(indxC);
    
    ChiWl(i,indxL)=1;
    ChiWl(i,indxR)=1;
    SChiWl        = SChiWl+ChiWl(i,:);
    
    SChi            = SChi+Chi(i,:);
    SshapeXcZc      = SshapeXcZc+shapeXcZc(i,:);
    SdxshapeXcZc    = SdxshapeXcZc+dxshapeXcZc(i,:);
    
    nu.z(i,:)=Chi(i,:);
    nu.x(i,:)=-dxshapeXcZc(i,:);
    nuz      =nuz+nu.z(i,:);
    nux      =nux+nu.x(i,:);
 end
Chi(i+1,:)=SChi;
ChiWl(i+1,:)=SChiWl;
shapeXcZc(i+1,:)=SshapeXcZc;
dxshapeXcZc(i+1,:)=SdxshapeXcZc;
nu.z(i+1,:)=nuz;
nu.x(i+1,:)=nux;
parS.shapeXcZc=shapeXcZc;parS.dxshapeXcZc=dxshapeXcZc;
parS.chi=Chi';parS.xShip=xShip;
parS.chiWl=ChiWl';
parS.nu=nu;

% eta=eta.*(1-chiship(:,end))+shippar.reform.shapeXcZc(end,:)';

% figure
% plot(x,SshapeXcZc,'b',x,parS.chi(:,end),'k--',x,SdxshapeXcZc,'r-.');
% title(['Mass (total)=',num2str(sum(Mass)),'; LWL (total)= ',num2str(sum(LWL))])
% plot_properties;

