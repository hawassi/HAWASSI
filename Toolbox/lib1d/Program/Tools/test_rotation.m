% dx=0.5;
% x=-20:dx:20;draft=5;
% z=zeros(size(x));
% indl=closest(x,-20);
% indr=closest(x,20);
% z(indl:indr)=-draft;
% dxz=zeros(size(x));
% dxz(indl)=-draft/dx;
% dxz(indr)=draft/dx;
% 
% theta=deg2rad(1);
% znow=z-theta.*(x+dxz.*z);
%  
% % figure;
% % plot(x,z,'r',x,znow,'b',x,dxz,'--k')
% % xlim([-30 30])


draft=5;
theta=deg2rad(10);
z0=[-draft:1:draft];
x0=[-20:1:20];
onez=ones(size(z0));
onex=ones(size(x0));

z=[fliplr(z0)        -draft.*onex z0            draft.*onex];
x=[x0(1)*onez x0          x0(end).*onez fliplr(x0)];
% figure;
% plot(x,z)

xnew=x.*cos(theta)-z.*sin(theta);
znew=x.*sin(theta)+z.*cos(theta);


figure;
fill(xnew,znew,'r')
%xlim([-25 25]);ylim([-10 10])
