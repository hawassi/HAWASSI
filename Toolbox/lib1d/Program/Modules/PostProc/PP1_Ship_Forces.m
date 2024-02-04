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
if strcmp(dynmodel(1:3),'HS1')
modelNonlin=1;
elseif strcmp(dynmodel(1:3),'HS2')
modelNonlin=2;    
elseif strcmp(dynmodel(1:3),'HS3')
modelNonlin=3; 
elseif strcmp(dynmodel(1:3),'HS4')
modelNonlin=4; 
end

indxi=closest(x,Xinterv(1)); 
indxf =closest(x,Xinterv(end));
indti=closest(time,tinterv(1)); 
indtf =closest(time,tinterv(2));
timeNow=time(indti:dplot:indtf);
Nt=length(timeNow);
g=9.81;
iter=1;Niter=Nt;tic;
set(jProgressBar,'Maximum',Niter, 'Value',iter);
%statusbarObj.setText(['estimating time remaining..']);
k=freqspace(x);


if any(shipsavevar.time(2:end)==shipsavevar.time(1))
IndEnd=find(shipsavevar.time(2:end)==shipsavevar.time(1),1)-1;
else
IndEnd=length(shipsavevar.time)-1; 
end
Nship=shippar.Nship;
velZ=zeros(length(timeNow),Nship);
velX=velZ;velTheta=velZ;
gammaZ=zeros(length(timeNow),Nship);
gammaX=gammaZ;gammaTheta=gammaZ;
FlinX=zeros(length(timeNow),Nship);
FlinZ=zeros(length(timeNow),Nship);
FlinThe=zeros(length(timeNow),Nship);


FnonlinX=zeros(length(timeNow),Nship);
FnonlinZ=FnonlinX;FnonlinTheta=FnonlinX;
dxi_KdiffZ=zeros(length(timeNow),Nship);
dxi_KdiffX=zeros(length(timeNow),Nship);
dxi_KdiffTheta=zeros(length(timeNow),Nship);

set(jProgressBar,'Maximum',Niter, 'Value',Niter/10);
%statusbarObj.setText(['estimating time remaining..']);

   
for ii=1:Nship
if ~strcmp(shippar.data(ii,2),'Fixed')
velZ(:,ii)=interp1(shipsavevar.time(1:IndEnd),shipsavevar.zVel(1:IndEnd,ii),timeNow);
velX(:,ii)=interp1(shipsavevar.time(1:IndEnd),shipsavevar.xVel(1:IndEnd,ii),timeNow);
velTheta(:,ii)=interp1(shipsavevar.time(1:IndEnd),shipsavevar.thetaVel(1:IndEnd,ii),timeNow);
aMvelZ=velZ(:,ii).*shippar.rad.Ma.z(ii);
aMvelX=velX(:,ii).*shippar.rad.Ma.x(ii);
else
aMvelZ=0;aMvelX=0;    
end
betaZ=interp1(shipsavevar.time(1:IndEnd),shipsavevar.betaz(1:IndEnd,ii),timeNow);
betaX=interp1(shipsavevar.time(1:IndEnd),shipsavevar.betax(1:IndEnd,ii),timeNow);
betaTheta=interp1(shipsavevar.time(1:IndEnd),shipsavevar.betatheta(1:IndEnd,ii),timeNow);
gammaTheta(:,ii)=betaTheta+aMvelZ.*velTheta(:,ii);
gammaZ(:,ii)=betaZ+aMvelZ.*velZ(:,ii);
gammaX(:,ii)=betaX+aMvelX.*velX(:,ii);


try
    BB=str2num(cell2mat(shippar.data(ii,3)));
    DD=str2num(cell2mat(shippar.data(ii,4)));
catch
    BB=cell2mat(shippar.data(ii,3));
    DD=cell2mat(shippar.data(ii,4));
end
factZ=9.811.*BB.*influx.Hs/2;
factX=9.811.*DD.*influx.Hs/2;
factThe=9.811.*BB.*DD.*influx.Hs/2;

FlinX(:,ii)= -gradient(gammaX(:,ii),timeNow(2)-timeNow(1))./factX;
FlinZ(:,ii)= -gradient(gammaZ(:,ii),timeNow(2)-timeNow(1))./factZ;
FlinThe(:,ii)=-gradient(gammaTheta(:,ii),timeNow(2)-timeNow(1))./factThe;

if isfield(shipsavevar,'dxi_Kdiffz')
dxi_KdiffZ(:,ii)=interp1(shipsavevar.time(1:IndEnd),shipsavevar.dxi_Kdiffz(1:IndEnd,ii),timeNow);
dxi_KdiffX(:,ii)=interp1(shipsavevar.time(1:IndEnd),shipsavevar.dxi_Kdiffx(1:IndEnd,ii),timeNow);
dxi_KdiffTheta(:,ii)=interp1(shipsavevar.time(1:IndEnd),shipsavevar.dxi_Kdifftheta(1:IndEnd,ii),timeNow);
FnonlinX(:,ii)=-dxi_KdiffX(:,ii)./factX;
FnonlinZ(:,ii)=-dxi_KdiffZ(:,ii)./factZ;
FnonlinTheta(:,ii)=-dxi_KdiffTheta(:,ii)./factThe;
end
end

 set(jProgressBar,'Maximum',Niter, 'Value',Niter/2);
 
    %%%%%%saving data 
    PPdata.PP1.ship_Force.FlinX=FlinX;
    PPdata.PP1.ship_Force.FlinZ=FlinZ;
    PPdata.PP1.ship_Force.FlinFlinThe=FlinThe;
    PPdata.PP1.ship_Force.FnonlinX=FnonlinX;
    PPdata.PP1.ship_Force.FnonlinZ=FnonlinZ;
    PPdata.PP1.ship_Force.FnonlinThe=FnonlinTheta;
    PPdata.PP1.ship_Force.time=timeNow;
    %%%%%%%%%%%%%%%

    hf6=figure('Name','PostProc: Forces', 'Position',[250,10,800,500]); 
    
    idCol=1;
    CR={'r','--g','-.b',':m','or','sg','db','+m'};
    Tp=influx.Tp;
    subplot(3,2,1)
    for ii=1:Nship
        if idCol==9, idCol=1; end;
        plot(timeNow./Tp,FlinX(:,ii),CR{idCol});
        hold on;
        idCol=idCol+1;
        minFlinX=min(FlinX(:,ii));
        maxFlinX=max(FlinX(:,ii));
        AFlinX=(maxFlinX-minFlinX)/2;
        leg{ii}=['Ship #',num2str(ii), ': (min,max,Av)=(',num2str(minFlinX,2),',',num2str(maxFlinX,2),',',num2str(AFlinX,2)];
    end
   
    title('Surge Force (dynamic)')
    xlabel('t/Tp') ;   ylabel('Fx/\rhogTA');
   % if Nship>1
        legend(leg{1:end},'Location','northwest');
   % end
    tinterv=GUIpp.PP1_Buoy_T/influx.Tp;
    if GUIpp.PPsetting_xlim==1
        xlim([GUIpp.PPsetting_xlimVal(1) GUIpp.PPsetting_xlimVal(2)])
    else
        xlim([tinterv(1) tinterv(2)])
    end
    
    if GUIpp.PPsetting_ylim==1
        ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
    end
    plot_properties;
    idCol=1;
     subplot(3,2,2)
    for ii=1:Nship
        if idCol==9, idCol=1; end;
        plot(timeNow./Tp,FnonlinX(:,ii),CR{idCol});
        hold on;
        idCol=idCol+1;
        minF=min(FnonlinX(:,ii));
        maxF=max(FnonlinX(:,ii));
        AF=(maxF-minF)/2;
        leg{ii}=['Ship #',num2str(ii), ': (min,max,Av)=(',num2str(minF,2),',',num2str(maxF,2),',',num2str(AF,2)];
    end
    title('Surge Force (nonlinear)')
    xlabel('t/Tp') ;    ylabel('Fx/\rhogTA');
   % if Nship>1
        legend(leg{1:end},'Location','northwest');
   % end
    
    if GUIpp.PPsetting_xlim==1
        xlim([GUIpp.PPsetting_xlimVal(1) GUIpp.PPsetting_xlimVal(2)])
    else
        xlim([tinterv(1) tinterv(2)])
    end
    
    if GUIpp.PPsetting_ylim==1
        ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
    end
    plot_properties;
    
     
    
    idCol=1;
     subplot(3,2,3)
    for ii=1:Nship
        if idCol==9, idCol=1; end
        plot(timeNow./Tp,FlinZ(:,ii),CR{idCol});
        hold on;
        idCol=idCol+1;
        minF=min(FlinZ(:,ii));
        maxF=max(FlinZ(:,ii));
        AF=(maxF-minF)/2;
        leg{ii}=['Ship #',num2str(ii), ': (min,max,Av)=(',num2str(minF,2),',',num2str(maxF,2),',',num2str(AF,2)];

    end
    title('Heave Force (dynamic)')
    xlabel('t/Tp') ;   ylabel('Fz/\rhogBA');
 %   if Nship>1
        legend(leg{1:end},'Location','northwest');
 %   end
    if GUIpp.PPsetting_xlim==1
        xlim([GUIpp.PPsetting_xlimVal(1) GUIpp.PPsetting_xlimVal(2)])
    else
        xlim([tinterv(1) tinterv(2)])
    end
    
    if GUIpp.PPsetting_ylim==1
        ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
    end
    plot_properties;
    idCol=1;
     subplot(3,2,4)
    for ii=1:Nship
        if idCol==9, idCol=1; end;
        plot(timeNow./Tp,FnonlinZ(:,ii),CR{idCol});
        hold on;
        idCol=idCol+1;
        minF=min(FnonlinZ(:,ii));
        maxF=max(FnonlinZ(:,ii));
        AF=(maxF-minF)/2;
        leg{ii}=['Ship #',num2str(ii), ': (min,max,Av)=(',num2str(minF,2),',',num2str(maxF,2),',',num2str(AF,2)];

    end
    title('Heave Force (nonlinear)')
    xlabel('t/Tp') ;   ylabel('Fz/\rhogBA');
 %   if Nship>1
        legend(leg{1:end},'Location','northwest');
 %   end
    if GUIpp.PPsetting_xlim==1
        xlim([GUIpp.PPsetting_xlimVal(1) GUIpp.PPsetting_xlimVal(2)])
    else
        xlim([tinterv(1) tinterv(2)])
    end
    
    if GUIpp.PPsetting_ylim==1
        ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
    end
    plot_properties;
    
   
    
    idCol=1;
     subplot(3,2,5)
    for ii=1:Nship
        if idCol==9, idCol=1; end;
        plot(timeNow./Tp,FlinThe(:,ii),CR{idCol});
        hold on;
        idCol=idCol+1;
        minF=min(FlinThe(:,ii));
        maxF=max(FlinThe(:,ii));
        AF=(maxF-minF)/2;
        leg{ii}=['Ship #',num2str(ii), ': (min,max,Av)=(',num2str(minF,2),',',num2str(maxF,2),',',num2str(AF,2)];

    end
    title('Pitch Force (dynamic)')
    xlabel('t/Tp') ;   ylabel('Ftheta/\rhogBTA');
    %if Nship>1
        legend(leg{1:end},'Location','northwest');
   % end
    if GUIpp.PPsetting_xlim==1
        xlim([GUIpp.PPsetting_xlimVal(1) GUIpp.PPsetting_xlimVal(2)])
    else
        xlim([tinterv(1) tinterv(2)])
    end
    
    if GUIpp.PPsetting_ylim==1
        ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
    end
    plot_properties;
    idCol=1;
     subplot(3,2,6)
    for ii=1:Nship
        if idCol==9, idCol=1; end;
        plot(timeNow./Tp,FnonlinTheta(:,ii),CR{idCol});
        hold on;
        idCol=idCol+1;
         minF=min(FnonlinTheta(:,ii));
        maxF=max(FnonlinTheta(:,ii));
        AF=(maxF-minF)/2;
        leg{ii}=['Ship #',num2str(ii), ': (min,max,Av)=(',num2str(minF,2),',',num2str(maxF,2),',',num2str(AF,2)];

    end
    title('Pitch Force (nonlinear)')
    xlabel('t/Tp') ;    ylabel('Ftheta/\rhogBTA');
    
 %   if Nship>1
        legend(leg{1:end},'Location','northwest');
 %   end
    
    if GUIpp.PPsetting_xlim==1
        xlim([GUIpp.PPsetting_xlimVal(1) GUIpp.PPsetting_xlimVal(2)])
    else
        xlim([tinterv(1) tinterv(2)])
    end
    
    if GUIpp.PPsetting_ylim==1
        ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
    end
    plot_properties;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     hf6b=figure('Name','PostProc: Total Forces and Spectra', 'Position',[250,10,800,500]); 
    
     BandFilt=GUIpp.Filter;
    if BandFilt==1
        LF=GUIpp.LFHF(1);      %Low Frequency
        HF=GUIpp.LFHF(2);      %High Frequency
    end
    
     
     subplot(3,2,1)
     idCol=1;
    for ii=1:Nship
        if idCol==9, idCol=1; end
        simul=FlinX(:,ii)+FnonlinX(:,ii);
        simul(isnan(simul))=0;
         if BandFilt==1
        simul=fun_bandpass(timeNow,simul,LF,HF);
         end
        plot(timeNow./Tp, simul,CR{idCol});
        hold on;
        idCol=idCol+1;
        minF=min( simul);
        maxF=max( simul);
        AF=(maxF-minF)/2;
        leg{ii}=['Ship #',num2str(ii), ': (min,max,Av)=(',num2str(minF,2),',',num2str(maxF,2),',',num2str(AF,2)];

    end
    title(' Surge Force (total)')
    xlabel('t/Tp') ;    ylabel('Fx/\rhogTA');
   % if Nship>1
        legend(leg{1:end},'Location','northwest');
  %  end
    
    if GUIpp.PPsetting_xlim==1
        xlim([GUIpp.PPsetting_xlimVal(1) GUIpp.PPsetting_xlimVal(2)])
    else
        xlim([tinterv(1) tinterv(2)])
    end
    
    if GUIpp.PPsetting_ylim==1
        ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
    end
    plot_properties;
    
      subplot(3,2,2)
     idCol=1;
     for ii=1:Nship
         if idCol==9, idCol=1; end
         simul=FlinX(:,ii)+FnonlinX(:,ii);
         simul(isnan(simul))=0;
          if BandFilt==1
         simul=bandpass(timeNow,simul,LF,HF);
          end
         [ww,spSIMUL] = spect(timeNow,simul,1,1);
         plot(ww,spSIMUL,CR{idCol});
         hold on;
         idCol=idCol+1;
         leg{ii}=['Ship #',num2str(ii)];
     end
    
    title('Surge Force (var. dens spect)')
    xlabel('\omega [rad/s]') ;    ylabel('Spect ');
   % if Nship>1
        legend(leg{1:end},'Location','northwest');
  %  end
    
%     if GUIpp.PPsetting_xlim==1
%         xlim([GUIpp.PPsetting_xlimVal(1) GUIpp.PPsetting_xlimVal(2)])
%     else
        xlim([ww(1) ww(floor(length(ww)/3))])
%     end
    
    if GUIpp.PPsetting_ylim==1
        ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
    end
    plot_properties;
    
    subplot(3,2,3)
    idCol=1;
    for ii=1:Nship
        if idCol==9, idCol=1; end;
         simul=FlinZ(:,ii)+FnonlinZ(:,ii);
         simul(isnan(simul))=0;
          if BandFilt==1
        simul=bandpass(timeNow,simul,LF,HF);
          end
        plot(timeNow./Tp,simul,CR{idCol});
        hold on;
        idCol=idCol+1;
         minF=min(simul);
        maxF=max(simul);
        AF=(maxF-minF)/2;
        leg{ii}=['Ship #',num2str(ii), ': (min,max,Av)=(',num2str(minF,2),',',num2str(maxF,2),',',num2str(AF,2)];

    end
    title('Heave Force (total)')
    xlabel('t/Tp') ;   ylabel('Fz/\rhogBA');
   % if Nship>1
        legend(leg{1:end},'Location','northwest');
   % end
    if GUIpp.PPsetting_xlim==1
        xlim([GUIpp.PPsetting_xlimVal(1) GUIpp.PPsetting_xlimVal(2)])
    else
        xlim([tinterv(1) tinterv(2)])
    end
    
    if GUIpp.PPsetting_ylim==1
        ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
    end
    plot_properties;
    
     subplot(3,2,4)
     idCol=1;
     for ii=1:Nship
         if idCol==9, idCol=1; end
         simul=FlinZ(:,ii)+FnonlinZ(:,ii);
         simul(isnan(simul))=0;
          if BandFilt==1
          simul=bandpass(timeNow,simul,LF,HF);
          end
         [ww,spSIMUL] = spect(timeNow,simul,1,1);
         plot(ww,spSIMUL,CR{idCol});
         hold on;
         idCol=idCol+1;
         leg{ii}=['Ship #',num2str(ii)];
     end
    
    title('Heave Force (var. dens spect)')
    xlabel('\omega [rad/s]') ;    ylabel('Spect ');
   % if Nship>1
        legend(leg{1:end},'Location','northwest');
  %  end
    
%     if GUIpp.PPsetting_xlim==1
%         xlim([GUIpp.PPsetting_xlimVal(1) GUIpp.PPsetting_xlimVal(2)])
%     else
        xlim([ww(1) ww(floor(length(ww)/3))])
%     end
    
    if GUIpp.PPsetting_ylim==1
        ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
    end
    plot_properties;
   
    idCol=1;
     subplot(3,2,5)
    for ii=1:Nship
        if idCol==9, idCol=1; end;
         simul=FlinThe(:,ii)+FnonlinTheta(:,ii);
         simul(isnan(simul))=0;
          if BandFilt==1
        simul=bandpass(timeNow,simul,LF,HF);
          end
        plot(timeNow./Tp,simul,CR{idCol});
        hold on;
        idCol=idCol+1;
        minF=min(simul);
        maxF=max(simul);
        AF=(maxF-minF)/2;
        leg{ii}=['Ship #',num2str(ii), ': (min,max,Av)=(',num2str(minF,2),',',num2str(maxF,2),',',num2str(AF,2)];

    end
    title('Pitch Force (total)')
    xlabel('t/Tp') ;    ylabel('Ftheta/\rhogBTA');
    
  %  if Nship>1
        legend(leg{1:end},'Location','northwest');
  %  end
    
    if GUIpp.PPsetting_xlim==1
        xlim([GUIpp.PPsetting_xlimVal(1) GUIpp.PPsetting_xlimVal(2)])
    else
        xlim([tinterv(1) tinterv(2)])
    end
    
    if GUIpp.PPsetting_ylim==1
        ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
    end
    plot_properties;

     subplot(3,2,6)
     idCol=1;
     for ii=1:Nship
         if idCol==9, idCol=1; end
         simul=FlinThe(:,ii)+FnonlinTheta(:,ii);
         simul(isnan(simul))=0;
          if BandFilt==1
          simul=bandpass(timeNow,simul,LF,HF);
          end
         [ww,spSIMUL] = spect(timeNow,simul,1,1);
         plot(ww,spSIMUL,CR{idCol});
         hold on;
         idCol=idCol+1;
         leg{ii}=['Ship #',num2str(ii)];
     end
    
    title('Pitch Force (var. dens spect)')
    xlabel('\omega [rad/s]') ;    ylabel('Spect ');
   % if Nship>1
        legend(leg{1:end},'Location','northwest');
  %  end
    
%     if GUIpp.PPsetting_xlim==1
%         xlim([GUIpp.PPsetting_xlimVal(1) GUIpp.PPsetting_xlimVal(2)])
%     else
        xlim([ww(1) ww(floor(length(ww)/3))])
%     end
    
    if GUIpp.PPsetting_ylim==1
        ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
    end
    plot_properties;
    set(jProgressBar,'Maximum',Niter, 'Value',Niter);
