function minR=funC_minradius(x,y)
xc=mean(x);
yc=mean(y);
minR=min(sqrt((x-xc).^2+(y-yc).^2));   
% figure(111)
% plot(x,y,'r',xc,yc,'or')