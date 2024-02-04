% mainABflat, version May 2011
% nsqrspec
% % % % Software developed by 
% % % % LabMath-Indonesia & University of Twente
%%
function [ww, nsqrspec, variance,gem]=nsqrspec(time,signal,pp)
    dt           =   time(2)-time(1);
    gem         =   mean(signal);
    variance    =   var(signal,1); %sum(signal.*signal)*dt/length(time);
    shat        =   abs(fft(signal));
    nsqrsp      =   sqrt(shat.*shat/sum(shat.*shat));
    ww          =   freqspace(time);
    nsqrspec    =   circshift(nsqrsp,floor((-pp/2)));
    nsqrspec    =   filter(ones(1,pp)/pp,1,nsqrspec);%circshift(nsqrspec,-floor(pp-1)));
    nsqrspec    =   (nsqrspec + fliplr(nsqrspec')')/2; 
   