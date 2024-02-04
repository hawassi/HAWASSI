function dat_pol=funC_CarttoPolCordinate(x,y,dat,r,theta)

[X,Y]     = meshgrid(x,y);    % meshgrid of cartesian data
[R,Theta] = meshgrid(r,theta);% meshgrid of polar data

X1=R.*cos(Theta);
Y1=R.*sin(Theta);
% figure;mesh(X,Y,dat);set(gcf,'renderer','zbuffer');view([0 90])
dat_pol=interp2(X,Y,dat,X1,Y1,'linear',0);