function [chi,chiPol]=funB_charcircle2D(x,y,xc,yc,r0,r1,theta,side,N)
dx=x(2)-x(1);dy=y(2)-y(1);
dr=sqrt(dx^2+dy^2);
r=0:dr:r1;
a=r0-N*(x(2)-x(1));
b=r0+N*(x(2)-x(1));

switch side
    case 'right'
%         char1D=funC_Heaviside(r-r0);
        char1D=funC_Heaviside(r-a).*sin((r-a)*pi/(b-a)/2).^2.*funC_Heaviside(b-r)+funC_Heaviside(r-b);
    case 'left'
%         char1D=funC_Heaviside(r0-r);
        char1D=1-(funC_Heaviside(r-a).*sin((r-a)*pi/(b-a)/2).^2.*funC_Heaviside(b-r)+funC_Heaviside(r-b));
end

chiPol=zeros(length(theta),length(r));
for j=1:length(theta);
    chiPol(j,:)=char1D;
end
chi=Fun_PoltoCartordinary(r,theta,chiPol,x,y,xc,yc,0);