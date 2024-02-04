function [param]=parambodysetup(dom,inputbody)

xyz0=inputbody.xyz0;
TT =inputbody.draft;
Nb  =inputbody.N;
Mass=zeros(Nb,1);
shapeN=0;chiN=0;
for i=1:Nb
[chiB]=shapeCharLib2d(dom,inputbody,i);   
 Mass(i)     = -trapz(dom.Y,trapz(dom.X,chiB.*(-TT(i)),2));
 shapeN=shapeN+chiB*(-TT(i)+xyz0(i,3));
 param.chiB(i).char=chiB;
 chiN=chiN+chiB;
end
param.chiB(Nb+1).char=chiN;
param.mass=Mass;
param.shape=shapeN;
[param.dxshape,param.dyshape]=gradient(param.shape,dom.dx,dom.dy);
hf=figure;
set(hf,'Renderer','zbuffer'); %due to graphics driver
subplot(2,1,1)
surf(dom.XX,dom.YY,chiN,'edgecolor','none')
%xlim([20 40])
subplot(2,1,2)
surf(dom.XX,dom.YY,shapeN,'edgecolor','none')
%xlim([20 40])

