function [omega,S]=analyse_cross_spectrum(time,x,y,domega)

% syntax:       [omega,S]=analyse_cross_spectrum(time,x,y,domega)
% example:      [omega,S]=analyse_cross_spectrum(time,x,y)
% description:  Function analyse_cross_spectrum estimates the cross
%               spectrum density of the signals x(time) and y(time), using
%               Welch method (half overlapping hanning windows)
%
% input:        - the time axis of signals x(time) and y(time) [s]
%               - the signals x and y
%               - optional, the desired domega stepp. If not specified
%               0.02rad/s will be used
%
% output        - omega [rad/s], the frequency axis
%               - S, the cross spectral density
%
% revisions:    Nr.     Date        By              Version
%                0    May 28 2010   JvdB               1.0 

s=size(x);
transpose=0;
if (s(1)>s(2))
    transpose=1;
    x=x';
end
s=size(y);
if (s(1)>s(2))
    y=y';
end
clear s;
if nargin < 4, domega= 0.02; end
dt=(time(end)-time(1))/(length(time)-1);
%uitdunning=floor(0.5/dt);
%if uitdunning>1,
%    disp(['uitdunning ' int2str(uitdunning)]);
%    t=t(1:uitdunning:end);
%    x=x(1:uitdunning:end);
%    dt=(t(end)-t(1))/(length(t)-1);
%end
m=pi/(dt*domega);
m=2*floor(m/2);
%disp([int2str(m) ' lags']);
[S,df]=crossSpectrumlags(x,y,dt,m);
S=S.';
f=(0:(length(S)-1))*df;
if transpose==1,
    f=f';
end
omega=f;

function [specx,df]=crossSpectrumlags(x,y,dt,m)
%
%
%    Input:   x     - Input signal
%             y     - Input signal 2
%             dt    - Step value
%             m     - Number of asked lags in spectrum
%                     auto-spectra (length m+1)
%
%    Output:  specx - Calculated spectrum
%             df    - Frequency step (rad/s)
%
%    Note:    2m may not be greater than n
%
%   December 1999 by R.F. Tap
%   Modified by J. van den Berg 2008, to compute crossSpectrum instead
%   of autoSpectrum and added comments
%

%
% Initialisations
%
n = length(x);
len = 2*m;
m1 = m+1;
specx(m1) = 0;
n1 = 0;
pi2 = 2*pi;
%
% Construeer halfoverlappende deelreeksen met lengte 2m en gemiddelde nul
%
md = floor(n/m)-1;
fac = dt/(pi*len);
for l=1:md
   xl = x(1+n1:len+n1);
   xl = xl-mean(xl);
   yl = y(1+n1:len+n1);
   yl = yl-mean(yl);
   n1 = n1+m;
   
   %
   % Data venster (Hanning window, dit gebeurt ook in Sally
   % de optionele Hanning smoothing in Sally is een extra post-hanning
   % (die eigenlijk zeer onwenselijk is, dus die zit hier niet in)
   %
   xl = 0.5*xl.*(1-cos(0:pi2/(len-1):pi2));
   yl = 0.5*yl.*(1-cos(0:pi2/(len-1):pi2));
   %
   % Fourier transformatie van deelreeksen
   %
   fx = fft(xl);
   fy = fft(yl);
   %
   % Periodogramschatter voor ieder deelinterval
   %
   specx = specx+fac*(fy(1:m1).*conj(fx(1:m1)));
end
%
% Bepaal spectrum door middeling van de periodogramschatters
%
specx = specx/md*8/3;
df = pi/(m*dt);