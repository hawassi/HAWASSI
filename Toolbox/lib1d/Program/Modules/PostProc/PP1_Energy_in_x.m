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

JavProgressBar;
[jbStop]=Java_stopbutton(statusbarObj);
jProgressBar.setStringPainted( true );



Xinterv=GUIpp.PP1_Prof_X';
tinterv=GUIpp.PP1_Buoy_T;
if GUIpp.PPsetting_tcoarse==1
    dplot= GUIpp.PPsetting_tcoarseVal;
else
    dplot=1;
end

if GUIpp.PPsetting_xcoarse==1
    dxplot= GUIpp.PPsetting_xcoarseVal;
else
    dxplot=1;
end
BandFilt=GUIpp.Filter;
if BandFilt==1
    LF=GUIpp.LFHF(1);      %Low Frequency
    HF=GUIpp.LFHF(2);     %High Frequency
end

indxi=closest(x,Xinterv(1)); 
indxf =closest(x,Xinterv(end));
indti=closest(time,tinterv(1)); 
indtf =closest(time,tinterv(2));
timeNow=time(indti:dplot:indtf);
Nt=length(timeNow);
ww=freqspace(time);
g=9.81;
iter=1;tic;
xNow=x(indxi:dxplot:indxf);
Niter=length(xNow);
set(jProgressBar,'Maximum',Niter, 'Value',iter);
statusbarObj.setText(['estimating time remaining..']);
Ep=zeros(size(xNow));
Ek=zeros(size(xNow));
IdFlat=0;
if min(bathy)==max(bathy)
 IdFlat=1;   
  Di=-bathy(1);
  kww=invOmExact(ww,-bathy(1));
  
  
  
end
dk=2*pi./(x(end)-x(1));
for i=1:length(xNow)
     IDstop=eventLoopStop(jbStop);
    if IDstop==1, break;end;
    
    indxX=closest(x,xNow(i));
    
    if BandFilt==1
       etaFilt=bandpass(time,eta(:,indxX),LF,HF);
       etaN=etaFilt(indti:dplot:indtf);
       uN=bandpass(time,u(:,indxX),LF,HF);
    else
       etaN=eta(indti:dplot:indtf,indxX);
       uN=u(:,indxX);
    end
    
    
    
   Ep(i) =trapz(timeNow,g*etaN.^2/2);  
   if IdFlat==0
   Di=-bathy(indxX);
   if Di<=0
      break; 
   end
   %%%preparing_for_interpolation;
   kmax=invOmExact(max(ww),Di);
   kmin=invOmExact(min(ww),Di);
   kk=linspace(kmin,kmax,floor((kmax-kmin)./dk));
   wwN=OmExact(kk,Di);
   %%% interpolation
   kww=interp1(wwN,kk,ww,'spline');
   %%%%%
   end
   Cpw=UpExact(kww,Di);
   K=0.5.*Ifft(Cpw.*fft(uN)).^2./g;
   Ek(i) =trapz(timeNow,K(indti:dplot:indtf));  
   
   if mod(iter,floor(0.1*Niter))==0 || iter==floor(0.01*Niter)
        set(jProgressBar,'Maximum',Niter, 'Value',iter);
        ETA=remain_time(iter,Niter);
        statusbarObj.setText(['time remaining=', num2str(ETA)]);
    end
    iter=iter+1;
    
end
 set(jProgressBar,'Maximum',Niter, 'Value',Niter);
statusbarObj.setText(['']);
jProgressBar.setVisible(0);
jbStop.setVisible(0)

figure;
plot(xNow,Ep+Ek,'r',xNow,Ep,'g',xNow,Ek,'b--');
xlabel('x [m]');
ylabel('Energy');
legend('Total','Potential','Kinetic')
plot_properties;
