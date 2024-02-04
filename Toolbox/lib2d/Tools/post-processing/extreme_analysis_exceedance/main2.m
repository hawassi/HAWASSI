clear
close all

%addpath ('i:\applications\scripts\maTLAB\');

%load 'D:\temp\ABCrest\140507 report_MARIN\data\data223002F_AB4U0_8.mat'
%load 'D:\Workspace\1. Codes\developer\Probability of Exceedance\140507 report_MARIN\data\data223002F_AB4U1.mat'
load 'D:\Workspace\1. Codes\developer\Probability of Exceedance\140507 report_MARIN\data\data223002F_AB4U0_8.mat'

%compare wave4LS @ X= 275.4
figure(2)
plot(data223002F.time,data223002F.simul(:,1),data223002F.time,data223002F.meas(:,1))
xlabel('time [s]')
ylabel('\eta [m]')
legend('AB','measured')
title('wave @ X=275.4')

%compare wave2_2L @ X= 312.4
figure(3)
plot(data223002F.time,data223002F.simul(:,2),data223002F.time,data223002F.meas(:,2))
xlabel('time [s]')
ylabel('\eta [m]')
legend('AB','measured')
title('wave @ X=312.4')

%compare wave6_4LS @ X= 642.6
figure(4)
plot(data223002F.time,data223002F.simul(:,3),data223002F.time,data223002F.meas(:,3))
xlabel('time [s]')
ylabel('\eta [m]')
legend('AB','measured')
title('wave @ X=642.6')

%compare wave13_LS @ X=  1593.4
figure(5)
plot(data223002F.time,data223002F.simul(:,4),data223002F.time,data223002F.meas(:,4))
xlabel('time [s]')
ylabel('\eta [m]')
legend('AB','measured')
title('wave @ X=1593.4')

for i=1:4
    %Analyze measured wave
    [waveheights_OB,T1_OB] = WaveExtremeAnalysis(data223002F.time',data223002F.meas(:,i)');

    %Determine probabilities of exceedance for measured wave heights
    [Probability_OB] = WaveExtremeProbability(waveheights_OB.CrestTrough,waveheights_OB.TroughCrest,waveheights_OB.Crest,waveheights_OB.Trough);

    %Analyze simulated wave
    [waveheights_AB,T1_AB] = WaveExtremeAnalysis(data223002F.time',data223002F.simul(:,i)');
    
    %Determine probabilities of exceedance for measured wave heights
    [Probability_AB] = WaveExtremeProbability(waveheights_AB.CrestTrough,waveheights_AB.TroughCrest,waveheights_AB.Crest,waveheights_AB.Trough);

    %Plot distributions on log scale
    figure(i)
    semilogy(waveheights_OB.Crest/waveheights_OB.SignificantHm0,Probability_OB.Crest,'bo')
    hold on
    semilogy(waveheights_AB.Crest/waveheights_AB.SignificantHm0,Probability_AB.Crest,'k*')
    hold off
    xlabel('\zeta_{\rm crest}/H_s [-]')
    ylabel('Probability of exceedance [-]')
    legend('measured [m]','AB [m]','Location','SouthWest')
    title(['wave probe at X= ' num2str(data223002F.waveprobes(i))])
end
