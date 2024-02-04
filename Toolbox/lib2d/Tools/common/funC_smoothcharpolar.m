function [chi,chiPol]=funC_smoothcharpolar(x,y,xc,yc,r0,r1,side)
dr=sqrt((x(2)-x(1))^2+(y(2)-y(1))^2);
r=linspace(0,r1,round(r1/dr));
dtheta=deg2rad(0.3);
theta=0:dtheta:pi;
switch side
    case 'right'
        char1D=heaviside(r-r0).*sin((r-r0)*pi/(r1-r0)/2).^2.*heaviside(r1-r)+heaviside(r-r1);
    case 'left'
        char1D=1-(heaviside(r-r0).*sin((r-r0)*pi/(r1-r0)/2).^2.*heaviside(r1-r)+heaviside(r-r1));
end

chiPol=zeros(length(theta),length(r));
for j=1:length(theta);
    chiPol(j,:)=char1D;
end
chi=funC_PoltoCartCordinate(r,theta,chiPol,x-xc,y-yc,0);