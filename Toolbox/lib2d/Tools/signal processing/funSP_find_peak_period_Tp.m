function Tp=funSP_find_peak_period_Tp(time,simul,Spsmooth)

ww= freqspace(time);
spSIMUL=abs(fft(simul));spSIMUL(1)=0;
if Spsmooth>1
    spSIMUL=smooth(spSIMUL,Spsmooth) ;
end
[~,sppeakindexPP]=max(spSIMUL);
nupeak     =ww(sppeakindexPP);
Tp=2*pi/nupeak;
end
