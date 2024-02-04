%A=load('AccessChannel_4_linear_7_T01_simul.mat');
%%

elev=InfluxSignal(2:end,floor(end/2));
time=InfluxSignal(2:end,1);
omsig=funC_freqspace(time);
% 
% %Fast Fourier Transform

varSig=var(elev);
% h00=4*sqrt(varSig);
% varSig=h00^2/16
S_hat=abs(fft(elev)).^2;
S_hat=S_hat(1:floor(end/2));
ww=omsig(1:floor(end/2));

var0=trapz(ww,S_hat);
S_hat=S_hat*varSig/var0;
var00=trapz(ww,S_hat);
h0=4*sqrt(var00)
VarDensSpect2=[abs(ww) S_hat];
save VarDensSpect2 VarDensSpect2
figure(33)
plot(ww,S_hat)
xlim([0 5])
xlabel('\omega')
