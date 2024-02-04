function JSsig=Jonswap(time,Tpeak,Hs,optionJS);%  30July 2011
% %  Jonswap spectrum and signal for environmental waves
% Jonswap spectrum is energy spectrum; for signal additional random phase
%%
% Tpeak = 1;
% optionJS = 'random';
% HS  = 3; 
% time=[0:0.05:2000*Tpeak-0.05];
dt = time(2)-time(1);
% Tp  = 10;
omp = 2*pi/Tpeak;
% dt  = Tp/20;
% 
% Tend = 40*Tp;
% NN      =2* floor(Tend/dt/2);
% 
% t   = [0:dt:(NN-1)*dt];
NN=2*floor(length(time)/2);
w   = freqspace(time,NN);
% figure(1)
% plot(w)

om  = w(1:NN/2);
switch optionJS
    case 'randsig'
        phase = rand(length(om),1);
%         ramp=zeros(NN,1);
        ramp=max(sin(pi/2*time/Tpeak).^2,sign(time-Tpeak))';
    case 'maxsig'
        phase = zeros(length(om),1);
        ramp=ones(NN,1);    
end; 
JSpec       = (om.^(-5)).*exp(-(5/4)*(omp./om).^4);
JShat      = sqrt(JSpec).*exp(2*pi*1i.*phase);
% figure(2)
% plot(om,abs(JShat));

JShat = [JShat; fliplr(conj(JShat(1:NN/2))')']*dt*pi/2;;%[JShat; 0;fliplr(conj(JShat(1:NN/2-1))')']*dt*pi/2;
% figure(3)
% plot(w,abs(JShat));

JSsig = ramp.*fftshift(ifft((JShat),'symmetric'));
varJSsig=var(JSsig);
hsquot = Hs/(4*sqrt(varJSsig));
JSsig = hsquot.* JSsig;
% figure(4) 
% plot(time,JSsig);
% figure(5)
% plot(w,abs(fft(JSsig)));
% 
% hscheck = 4*sqrt(var(JSsig))
end

%%