function bathyslope = Bathy_slope(x,parslope)
% parslope = [D_left,D_right,steepness,Xstartslope];

D_left      = parslope(1);              %depth_left 
D_right     = parslope(2);              %depth_right 
sf          = parslope(3);              % steepness fraction 
Xss         = parslope(4);              % start slope 
indXss      = closest(x,Xss);
length_s    = abs(D_right -D_left)/sf;    % length of slope
Xse         = Xss + length_s;
indXse      = closest(x,Xse);

%  bathyslope = ( D_left*ones(1,length(x)) + ...
%             -sign(D_left-D_right)*sf*(x-Xss).*Heaviside(x-Xss) ).*Heaviside(Xse-x)+ ...
%             + D_right.*Heaviside(x-Xse);

bathyslope =  D_left*ones(1,length(x));
bathyslope(indXss:indXse)=D_left-sign(D_left-D_right)*sf*(x(indXss:indXse)-Xss);
bathyslope(indXse+1:end)=D_right;
bathyslope=-bathyslope;
% figure;
% plot(x,bathyslope)
end
        

