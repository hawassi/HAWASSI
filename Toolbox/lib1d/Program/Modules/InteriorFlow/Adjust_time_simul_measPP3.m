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

function measAdj=Adjust_time_simul_measPP3(simul,meas,Zmeas)
measAdj=zeros(length(meas(:,1)),length(Zmeas));
for i=1:length(Zmeas)
    [Xcorel,lags]=xcorr(simul(:,i),meas(:,i),'coeff');
    OptCorr1=max(Xcorel);
    indMaxCorel=closest(Xcorel,OptCorr1);
    OptCorr2=lags(indMaxCorel);
    measAdj(:,i) = circshift(meas(:,i),OptCorr2);
end