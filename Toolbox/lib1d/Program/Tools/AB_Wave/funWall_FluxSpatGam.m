function [gamX_hat,gamX_skew_hat]=funWall_FluxSpatGam(par)
x=par.x;
k=par.k;
dx=par.dx;g=par.g;
Xwall=par.wall.position;
depth=-par.bathy(closest(x,Xwall));
Inarea=par.wall.gamX_Check;
gamX_hat        = zeros(length(k),1);
fact=1./dx;
if strcmp(Inarea, 'Area')
   
    vgx         = Ifft(UgExact(k,depth));
    dista       = x(closest(x,Xwall))-x(1);
    vgx         = circshift(vgx,floor(dista/dx));
    gamX_hat        = fft(vgx)*fact;
    gamX_hat        = gamX_hat*sqrt(g*depth)/gamX_hat(1)/dx;
    gamX_skew_hat   = 1i*OmExact(k,depth).*gamX_hat; % for uni-dir influxing 
    %gamX_skew_hat   = fft(smooth(Ifft(gamX_skew_hat),5));% 3*fft(gradient(Ifft(gamX_hat),dx));%
elseif strcmp(Inarea, 'Point')
    dista       = x(closest(x,Xwall))-x(1);
    gamX_hat      = exp(-1i*k*dista)*fact;
    gamX_skew_hat = 1i*k.*gamX_hat;
    gamX_skew_hat = fft(smooth(x,Ifft(gamX_skew_hat),0.01,'loess'));
    % Modify input signals...  
end

if strcmpi(par.wall.wallInfdir,'Uni+')
  gamX_skew_hat=-gamX_skew_hat;  
end  
%    figure;
%    plot(x,Ifft(gamX_skew_hat),'r',x,Ifft(gamX_hat),'b')
end