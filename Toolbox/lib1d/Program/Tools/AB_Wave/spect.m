% mainABflat, version May 2011
% nsqrspec
% % % % Software developed by
% % % % LabMath-Indonesia & University of Twente
%%
function [ww, spec, variance,gem]=spect(time,signal,pp,EspId)
dt           =   time(2)-time(1);
gem          =   mean(signal);
signal       =   signal-gem;
    
if EspId ==1
    shat        =   (abs(fft(signal))).^2;
else
    shat        =   abs(fft(signal));
end
variance     =   var(signal,1); %sum(signal.*signal)*dt/length(time);

ww      =   freqspace(time);
if pp==1
spec    =shat;    
else
% spec    =   circshift(shat,floor((-pp/2)));
% spec    =   filter(ones(1,pp)/pp,1,spec);%circshift(nsqrspec,-floor(pp-1)));
% spec    =   (spec + fliplr(spec')')/2;
spec=smooth(shat,pp);
end
ww=ww(1:floor(end/2));

if EspId==1
 var0=trapz(ww,spec(1:floor(end/2)));  
 spec=spec(1:floor(end/2))*variance/var0; %%only_half)side
%  [variance/var0 1/length(signal)];
%    varnow=trapz(ww,spec)
else
 spec=spec(1:floor(end/2));
end
