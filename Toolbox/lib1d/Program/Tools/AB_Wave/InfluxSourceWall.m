%% This function constructs source to be added to ode-solver

function [gamX,gamX_skew] = InfluxSourceWall(x,k,g,Xwall,Xleft,depth,Inarea,Om,Ug)
% influx.position     =  0;               % [m] (should be inside computation interval, not interfering with damping zones)
% influx.area         = 'GVArea';         % = 'GVArea' uses group-velocity area centered at .position with specified influx signal;
%                                         % = 'Point' uses delta-function influx  at .position with modified signal
%                                         % = 'User_influx'
% influx.direction       = 'Uni+';           % = 'Uni+' for unidrectional to right
%                                         % = 'Uni-' for unidrectional to left
%                                         % = 'Bisymm' for bidirectional symmetric both sides
% influx.nonlin_adjust  =  0;             % length (in # of  peak-wavelength) for adjustment-zone
% X_influx     = influx.position;


%%%%%%%%%%%%%%%%
%% Influx-spatial functions
%     Xinflux     = influx.position;
dx = x(2)-x(1);%(Xright-Xleft)/length(x);
gamX        = zeros(length(k),1);
fact        = 1/dx;%length(k)/(Xright-Xleft); %=1/dx

if strcmp(Inarea, 'GVArea')
    
    %         ug = Ug(:,depth);
    vgx         = Ifft(Ug(k,depth));
    dista       = Xwall - Xleft;%x(closest(x,Xinflux))-x(closest(x,Xleft));
    vgx         = circshift(vgx,floor(dista/dx));
    gamX        = fft(vgx)*fact;
    gamX        = gamX*sqrt(g*depth)/gamX(1)/dx;
    gamX_skew   = 1i*Om(k,depth).*gamX; % for uni-dir influxing
    gamX_skew   = fft(smooth(Ifft(gamX_skew),0.005,'loess'));
%     figure 
%     plot(x,Ifft(gamX),'r',x,Ifft(gamX_skew))
    
elseif strcmp(Inarea, 'Point')
    dista     = x(closest(x,Xwall))-x(closest(x,Xleft));
    gamX      = exp(-1i*k*dista)*fact;
    gamX_skew = 1i*Om(k,depth).*gamX;
    gamX_skew = fft(smooth(x,Ifft(gamX_skew),0.01,'loess'));
    % Modify input signals...  
end;


