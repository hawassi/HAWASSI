function [xx1,xx2,dx,Nx]=fun_obtainSimulSpaceShipInterval(xS,Ls,nS,xl,xr)
dx=Ls/(nS-2);
Sl=xS-dx*(nS-1)/2;
Sr=xS+dx*(nS-1)/2;
Nl=floor((Sl-xl)/dx+1);
Nr=floor((xr-Sr)/dx+1);
x1=Sl-Nl*dx;
x2=Sr+Nr*dx;
xx=x1:dx:x2;

fact=factor(length(xx));

idAdd=1;
while fact(end)>13

Nl=floor((Sl-xl)/dx+1);
Nr=floor((xr-Sr)/dx+1)+idAdd;
x1=Sl-Nl*dx;
x2=Sr+Nr*dx;
xx=x1:dx:x2;
idAdd=idAdd+1;
fact=factor(length(xx));
end
xx1=xx(1);xx2=xx(end);

Nx=length(xx);
ship=zeros(size(xx));
indl1=closest(xx,Sl);
indl2=closest(xx,Sr);
ship(indl1:indl2)=1;
% figure;
% plot(xx,ship)
