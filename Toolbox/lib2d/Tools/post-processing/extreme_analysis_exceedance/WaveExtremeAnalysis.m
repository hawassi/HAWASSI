function [Waveheights,T1] = WaveExtremeAnalysis(time,zeta)

% usage : [Waveheights] = WaveExtremeAnalysis(zeta);

% input:
% zeta = array with wave elevation [m]
% time = array with time [s]
%
% output:
% Waveheights.CrestTrough = waveheights crest to trough
% Waveheights.TroughCrest = waveheights crest to trough
% Waveheights.Trough = trough depths
% Waveheights.Crest = crest heights
% Waveheights.SignificantHm0 = Hs
% T1

WAVE_CL=zeta;

%compute wave difference to find extremes
dWAVE_CL_for=0*WAVE_CL;
dWAVE_CL_for(1:end-1)=WAVE_CL(2:end)-WAVE_CL(1:end-1);
dWAVE_CL_aft=0*WAVE_CL;
dWAVE_CL_aft(2:end)=WAVE_CL(2:end)-WAVE_CL(1:end-1);

%check where wave difference changes sign
FF=dWAVE_CL_aft.*dWAVE_CL_for;
index=find(FF <= 0);

%wave extremes
extremes=WAVE_CL(index);

%check where extremes change sign (zero-crossing)
S=sign(extremes);
Schange=S(1:end-1).*S(2:end);
I=find(Schange < 0);

%loop over extremes to determine max crest, max trough and max wave height
imin1=index(1);
k1=0;
k2=0;
for k=1:length(I)-1
   imax1=I(k);
   imin2=imax1+1;
   imax2=I(k+1);
   E1=extremes(imin1:imax1);
   E2=extremes(imin2:imax2);
   if (sum(E1) > 0)
       k1=k1+1;
       maximum=max(E1);
       minimum=min(E2);
       Hmax_crest_trough(k1)=maximum-minimum;
       Crest(k1)=maximum;
       Trough(k1)=minimum;
   else
       k2=k2+1;
       maximum=max(E2);
       minimum=min(E1);
       Hmax_trough_crest(k2)=maximum-minimum;
   end
   imin1=imin2;
end

%sort values
Hmax_crest_trough=sort(Hmax_crest_trough);
Hmax_trough_crest=sort(Hmax_trough_crest);
Crest=sort(Crest);
Trough=sort(-Trough);
Trough=-Trough;

Waveheights.CrestTrough=Hmax_crest_trough;
Waveheights.TroughCrest=Hmax_trough_crest;
Waveheights.Crest=Crest;
Waveheights.Trough=Trough;

%spectral moment
m=std(WAVE_CL)^2;
Waveheights.SignificantHm0=4*sqrt(m);
%!
[om,S]=analyse_auto_spectrum(time,WAVE_CL,0.01);
%!
T1=2*pi*sum(S)/sum(S.*om');

return
end
