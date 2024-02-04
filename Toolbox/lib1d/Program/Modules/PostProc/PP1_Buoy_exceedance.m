% !--------------------------------------------------------------------------------------
% !
% !    Copyright (C) 2024 - LabMath-Indonesia
% !
% !    This program is free software: you can redistribute it and/or modify
% !    it under the terms of the GNU General Public License as published by
% !    the Free Software Foundation, either version 3 of the License, or
% !    (at your option) any later version.
% !
% !    This program is distributed in the hope that it will be useful,
% !    but WITHOUT ANY WARRANTY; without even the implied warranty of
% !    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% !    GNU General Public License for more details.
% !
% !    You should have received a copy of the GNU General Public License
% !    along with this program.  If not, see <http://www.gnu.org/licenses/>.
% !
% !   Contributors list:
% !   - R. Kurnia
% !   - E. van Groesen
% !--------------------------------------------------------------------------------------



tinterv=GUIpp.PP1_Buoy_T;
ti=tinterv(1); tf=tinterv(2);
idt1=closest(time,ti);
idt2=closest(time,tf);

if GUIpp.PPsetting_tcoarse==1
    dplot= GUIpp.PPsetting_tcoarseVal;
else
    dplot=1;
end

Xbuoy   =GUIpp.PP1_Buoy_X';
buoyind = closest(x,Xbuoy);
etabuoy = eta(:,buoyind);

BandFilt=GUIpp.Filter;
if BandFilt==1
    LF=GUIpp.LFHF(1);      %Low Frequency
    HF=GUIpp.LFHF(2);     %High Frequency
end
idCol=1;
CR1={'r','g','b','m'};
CR2={'--r','--g','--b','--m'};
leg=[];
timeN=time(idt1:dplot:idt2);
hfEx=figure('Name','PostProc', 'Position',[250,10,800,500]);
 idLeg=1;
for I=1:length(Xbuoy)
    if mod(idCol,5)==0, idCol=1; end
    if BandFilt==1
        etaN=bandpass(timeN,etabuoy(idt1:dplot:idt2,I),LF,HF);
    else
        etaN=etabuoy(idt1:dplot:idt2,I);
    end
     %Analyze simulated wave
    
    [waveheights_AB_eta,T1_AB_eta] = WaveExtremeAnalysis(timeN,etaN');
    %Determine probabilities of exceedance for measured wave heights
    [Probability_AB_eta] = WaveExtremeProbability(waveheights_AB_eta.CrestTrough,waveheights_AB_eta.TroughCrest,waveheights_AB_eta.Crest,waveheights_AB_eta.Trough);
    
    [waveheights_AB_mineta,T1_AB_mineta] = WaveExtremeAnalysis(timeN,-etaN');
    %Determine probabilities of exceedance for measured wave heights
    [Probability_AB_mineta] = WaveExtremeProbability(waveheights_AB_mineta.CrestTrough,waveheights_AB_mineta.TroughCrest,waveheights_AB_mineta.Crest,waveheights_AB_mineta.Trough);
   
    
    %Plot distributions on log scale
    Hs_eta=waveheights_AB_eta.SignificantHm0;
%     Hs_mineta=waveheights_AB_mineta.SignificantHm0;

    semilogy(waveheights_AB_eta.Crest/Hs_eta,Probability_AB_eta.Crest,CR1{idCol},...
        waveheights_AB_mineta.Crest/Hs_eta,Probability_AB_mineta.Crest,CR2{idCol});
        PPdata.PP1.exceed(I).probCresty=Probability_AB_eta.Crest';
        PPdata.PP1.exceed(I).probCrestx=(waveheights_AB_eta.Crest/Hs_eta)';
        PPdata.PP1.exceed(I).probTroughy=Probability_AB_mineta.Crest';
        PPdata.PP1.exceed(I).probTroughx=(waveheights_AB_mineta.Crest/Hs_eta)';

    hold on;
    idCol=idCol+1;
    leg{idLeg}=[' Cr. exc. x=',num2str(Xbuoy(I)),'[m]'];
    leg{idLeg+1}=[' Tr. exc. x=',num2str(Xbuoy(I)),'[m]'];
    idLeg=idLeg+2; 
end


xlabel('|\eta|/H_{s}ref. [-]')
ylabel('Probability of exceedance [-]')
legend(leg{1:end},'Location','northwest');
plot_properties;
if GUIpp.save_figure==1
    saveas(hfEx,[sf_savename,'PP1_Buoy_exceedance','.fig'])
end