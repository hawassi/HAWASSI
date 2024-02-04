function [Probability] = WaveExtremeProbability(Hct,Htc,Cr,Tr)

% usage : [Probability] = WaveExtremeProbability(Hct,Htc,Cr,Tr)

% input:
% Hct = waveheights crest to trough
% Htc = waveheights trough to crest
% Cr = crest heights
% Tr = trough depths
% 
% output:
% Probability.CrestTrough = probability of exceedance waveheights crest to trough
% Probability.TroughCrest = probability of exceedance waveheights crest to trough
% Probability.Trough = probability of exceedance trough depths
% Probability.Crest = probability of exceedance crest heights

%sort values
Hmax_crest_trough=sort(Hct);
Hmax_trough_crest=sort(Htc);
Crest=sort(Cr);
assignin('base','Crest',Crest);
Trough=sort(-Tr);
Trough=-Trough;

I1=length(Hmax_crest_trough):-1:1;
I1=I1/length(Hmax_crest_trough);

I2=length(Hmax_trough_crest):-1:1;
I2=I2/length(Hmax_trough_crest);

I3=length(Crest):-1:1;
assignin('base','I3b',I3);
I3=I3/length(Crest);
assignin('base','I3f',I3);

I4=length(Trough):-1:1;
I4=I4/length(Trough);

Probability.CrestTrough=I1;
Probability.TroughCrest=I2;
Probability.Crest=I3;
Probability.Trough=I4;

return
end
