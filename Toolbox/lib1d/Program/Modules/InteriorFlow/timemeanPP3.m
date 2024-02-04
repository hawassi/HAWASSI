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

function [TMPhysVal_simul,TMPhysVal_meas]=timemeanPP3(simul,meas,Zmeas,D,Tsimul)
TMPhysVal_simul=zeros(length(Zmeas),1);
TMPhysVal_meas=zeros(length(Zmeas),1);
indt1=closest(Tsimul,Tsimul(1));
indt2=closest(Tsimul,Tsimul(end));

g=9.81;
for i=1:length(Zmeas)
    simulorg=simul(indt1:indt2,i);
    simulorg(isnan(simulorg))=0;
    [Xcorel,lags]=xcorr(simulorg,meas(indt1:indt2,i),'coeff');
    OptCorr1=max(Xcorel);
    indMaxCorel=closest(Xcorel,OptCorr1);
    OptCorr2=lags(indMaxCorel);
    measOpt = circshift(meas(indt1:indt2,i),OptCorr2);
    
    TMPhysVal_simul(i,1)=mean(simulorg)./sqrt(g*D);
    TMPhysVal_meas(i,1) =mean(measOpt)./sqrt(g*D);
end