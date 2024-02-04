% function SeaINxt =GenSeaIn(x,timesig,ptsignal,Omd)
% output is a function(x,t) to be used on (part of x) as influx source
% create omega-space, and amplitude spectrum
t=timesig;
ww = freqspace(t);
ampli = abs(fft(insig));
timephase = @(t)(exp(-1i*ww*t));
wwhalf = ww(2:floor(length(ww)/2));
%%
figure
plot(ww,ampli)
%%
% define the wavenumbers of x-space
k = freqspace(x);
dk = k(3)-k(1);
%%
% find wavenumbers corresponding to ww

kk = interp1(Omd(k),k,wwhalf);  %% these are the length of wavenumbers that correspond to ww
ind = isnan(kk);
kk(ind)=max(kk);
%%
figure
plot(kk,Omd(kk),'r',k,Omd(k)','b')
%%
% @(x,m) (    exp(1i*(kk(m)*cos(theta(m))*x+alp(m)))   )
% simpler
@(x,m) (   ampli(m)* exp(1i*(kk(m)))   )

