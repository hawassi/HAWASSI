%% This function constructs source to be added to ode-solver
function [INsig,gamX_hat,INsig_skew,gamX_skew_hat,insig,insig_skew] = ...
InfluxSource(par,time,infl_sig,infl_sigskew,Infmodel,Om,Up,Ug,om2k,GVom, lambda_p,omAdd,ugAdd)

x=par.x; k=par.k; g=par.g; depth=par.depth;
Xinflux=par.Xinflux; Xleft=par.Xleft;

sigin      = infl_sig;
sigin_skew  = infl_sigskew;

dx = x(2)-x(1);
gamX_hat        = zeros(length(k),1);gamX_skew_hat=gamX_hat;
fact        = 1/dx;%length(k)/(Xright-Xleft); %=1/dx

if strcmp(Infmodel,'GVArea')
    vgx         = Ifft(Ug(k,depth,ugAdd));
    dista       = Xinflux - Xleft;%x(closest(x,Xinflux))-x(closest(x,Xleft));
    vgx         = circshift(vgx,round(dista/dx));
    gamX_hat        = fft(vgx)*fact;
    gamX_hat        = gamX_hat*sqrt(g*depth)/gamX_hat(1)/dx;
    gamX_skew_hat   = 1i*Om(k,depth,omAdd).*gamX_hat; % for uni-dir influxing 1i*k.*gamX_hat;%
    gamX_skew_hat   = fft(smooth(Ifft(gamX_skew_hat),3));
    
elseif strcmp(Infmodel,'GVAreaShort')
    lsup=0.4.*lambda_p;Lb=2.*lsup.^2;
    gaussInf=zeros(length(x),1);
    ind1=closest(x,Xinflux-2.*lambda_p);
    ind2=closest(x,Xinflux+2.*lambda_p);
    gaussInf(ind1:ind2)=exp(-(x(ind1:ind2)-Xinflux).^2/Lb);  %gaussian function
    
    vgx       = Ifft(Ug(k,depth,ugAdd));
    dista     = Xinflux - Xleft;%x(closest(x,Xinflux))-x(closest(x,Xleft));
    vgx       = circshift(vgx,round(dista/dx)).*gaussInf;
    gamX_hat      = fft(vgx)*fact;
    gamX_hat      = gamX_hat*sqrt(g*depth)/gamX_hat(1)/dx;
    gamX_skew_hat = 1i*Om(k,depth,omAdd).*gamX_hat; % for uni-dir influxing 1i*k.*gamX_hat;%;%
    gamX_skew_hat = fft(smooth(Ifft(gamX_skew_hat),3));
    
    gamXom     = interp1(k,abs(gamX_hat),om2k,'spline');
    sigin      = Ifft(fft(sigin).*GVom'./gamXom');
    sigin_skew = cumtrapz(time,sigin);%Ifft(Up(om2k,depth)'.*fft(cumtrapz(time,sigin))); %integration of influx
elseif strcmp(Infmodel,'Point')
    dista    = x(closest(x,Xinflux))-x(closest(x,Xleft));
    gamX_hat     = exp(-1i*k*dista)*fact;
    gamX_skew_hat = 1i*Om(k,depth,omAdd).*gamX_hat;%1i*k.*gamX_hat;%
    gamX_skew_hat = fft(smooth(x,Ifft(gamX_skew_hat),0.01,'loess'));
end

%% Input signals
INsig       = spline(time,sigin);
INsig_skew  = spline(time,sigin_skew);
insig       = sigin;    %updated
insig_skew  = sigin_skew; 


