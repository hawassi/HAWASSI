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

Xinterv=GUIpp.PP1_Prof_X';
tinterv=GUIpp.PP1_Buoy_T;
if GUIpp.PPsetting_tcoarse==1
    dplot= GUIpp.PPsetting_tcoarseVal;
else
    dplot=1;
end

indxi=closest(x,Xinterv(1)); 
indxf =closest(x,Xinterv(end));
indti=closest(time,tinterv(1)); 
indtf =closest(time,tinterv(2));
timeNow=time(indti:dplot:indtf);
if ~isempty(GUIpp.PP1_Buoy_Spectrum_Smooth)
Spsmooth=GUIpp.PP1_Buoy_Spectrum_Smooth;
else
Spsmooth=1;   
end

HSig     =zeros(length(buoyind),1);
Tm01     =zeros(length(buoyind),1);
lambda_m01=zeros(length(buoyind),1);
Tpeak    =zeros(length(buoyind),1);
lambda_p =zeros(length(buoyind),1);
Vari     =zeros(length(buoyind),1);
Kurt     =zeros(length(buoyind),1);
Sk       =zeros(length(buoyind),1);
Asym     =zeros(length(buoyind),1);
% VarDensSpect=zeros(length(buoyind),1);

EspId=0;
for i=1:length(buoyind)
    simul=eta(indti:indtf,buoyind(i));
    %[ww,spSIMUL,varsimul] = spect(timeNow,simul,Spsmooth,EspId);
    ww= freqspace(time);  varsimul= var(simul,1);
    spSIMUL=abs(fft(simul));spSIMUL(1)=0;
    if Spsmooth>1
     spSIMUL=smooth(spSIMUL,Spsmooth) ;  
    end
    [~,sppeakindexPP]=max(spSIMUL);
    HSig(i)    =4*sqrt(varsimul);
    nupeak     =ww(sppeakindexPP); depth=-bathy(buoyind(i));
    kpeak      =invOmExact(nupeak,depth);
    Tm01(i)    =funSP_find_mean_period_Tm01(time(indti:indtf),simul);
    lambda_m01(i)=2*pi/invOmExact(2*pi/Tm01(i),depth);
    Tpeak(i)   = 2*pi/nupeak;
    lambda_p(i)=2*pi/kpeak;
    Vari(i)    =varsimul;
    Kurt(i)    =kurtosis(simul);
    MeanSim    =mean(simul);
    HilbSim    =imag(hilbert(simul-MeanSim));
    Sk(i)      =mean((simul-MeanSim).^3)/((mean((simul-MeanSim).^2))^(3/2));
    Asym(i)    =mean(HilbSim.^3)/((mean((simul-MeanSim).^2))^(3/2));
%     if EspId~=1
%         [ww,spSIMULEn] = spect(timeNow,simul,Spsmooth,1);
%         VarDensSpect(i)=trapz(spSIMULEn)*(ww(2)-ww(1));
%     else
%         VarDensSpect(i)=trapz(spSIMUL)*(ww(2)-ww(1));
%     end
end


%f = figure('Position',[440 500 461 146]);
f = figure('units','normalized','Position',[0.65 0.5 0.25 0.13],'name','Quantitative information.');

d = [HSig';Tm01';lambda_m01';Tpeak';lambda_p';Vari';Kurt';Sk';Asym'];
NBuoy=length(buoyind);
%        Create the column and row names in cell arrays
for i=1:NBuoy
    cnames{i} = ['x=',num2str(x(buoyind(i)),3)];
end

rnames = {'Hs','Tm01','lambda_m01','Tp','lambda_p','Var','Kurt','Sk','As'};
%      Create the uitable
t = uitable(f,'units','normalized','Position',[0 0 1 1],'Data',d,...
    'ColumnName',cnames,...
    'RowName',rnames,'ColumnWidth','auto');
title('Table. Quantity information')
if GUIpp.save_figure==1
    saveas(f,[sf_savename,'PP1_Buoy_quantity','.fig'])
end
