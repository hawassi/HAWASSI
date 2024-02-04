function [dat_cart,Theta1]=funC_PoltoCartCordinate(r,theta,dat,x,y,opt)

[X,Y]     = meshgrid(x,y);    % meshgrid of cartesian data
[R,Theta] = meshgrid(r,theta);% meshgrid of polar data
dtheta    =theta(2)-theta(1);
R1        =sqrt(X.^2+Y.^2);
Theta1    =atan2(Y,X);
% ind       =floor(length(y)/2);
[~,ind]   =min(Theta1(:,1));
if opt==0;
    Theta1(1:ind,:)=Theta1(1:ind,:)+(2*pi)*ones(length(1:ind),size(Theta1,2));
    % %range 0-2*pi
else
    Theta1(ind+1:end,1:floor(length(x)/2))=Theta1(ind+1:end,1:floor(length(x)/2))-2*pi*ones(length(ind+1:size(Theta1,1)),length(1:floor(length(x)/2)));
    % range peter : -270 degree till 90 degree
end

dat_cart=interp2(R,Theta,dat,R1,Theta1,'linear',0);

