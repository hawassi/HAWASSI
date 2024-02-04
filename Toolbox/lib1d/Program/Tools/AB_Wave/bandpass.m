function [Y] = bandpass(t,X,wfiltlo,wfiltup,plotswitch);
% BANDPASS.M
% 
% Band pass filter: Y = bandpass(t,X,wfilt1o,wfiltup)
%
% Input:  time base 't'
%         signal 'X'
%         filter frequency 'wfiltlo' [rad/s]: lower boundary (high pass filter)
%         filter frequency 'wfiltup' [rad/s]: upper boundary (low pass filter)
%
% Output: filtered signal 'Y'
%
%
% 10-jul-2006: modified to check on consistent length of t & X
%
% 11-sep-2006: modified to avoid spurious oscillations at ends of timetrace

t = t(:); X = X(:);

if nargin < 5,
    plotswitch = 0;
end
%
if length(t) ~= length(X), error('Vectors t & X must be the same lengths'), end
%
t = t-t(1);
%
N    = length(X);
dt   = (t(N) - t(1))/(N-1);
wNyq = pi/dt;	% Nyquist frequency = 1/2 sampling frequency
if wfiltup > wNyq,
   wfilt = wNyq;
   disp('warning: lowpass filter frequency > Nyquist frequency')
   disp('         lowpass filter frequency set to Nyquist frequency')
end
if wfiltlo >= wfiltup,
   wfiltlo = 0;
   disp('warning: highass filter frequency > lowpass filter frequency')
   disp('         highpass filter frequency set to zero')
end
% houtje-touwtje
tt=0:dt:(3*N-1)*dt;
N    = length(tt);
XX=[-X(end:-1:1)+2*X(1); X; -X(end:-1:1)+2*X(end)];
FX = fft(XX);
FY = FX;
dom= 2*pi/tt(N);
om = 0:dom:(N-1)*dom;
y1 = find(om > wfiltup);
y2 = find(om < wfiltlo);
FY(y1) = 0;
FY(y2) = 0;
N2 = rnd(N/2,0);
i  = 0:N2;
FY(N-i) = conj(FY(i+2));
YY  = real(ifft(FY));
Y=YY(end/3+1:end/3*2);
%
if plotswitch == 1,
    plot(om,abs(FX),om,abs(FY),'r')
    pause
    close
end
return
% official method
w1    = wfilt/wNyq;
B     = fir1(n,w1);
Y     = filtfilt(B,1,X);
return