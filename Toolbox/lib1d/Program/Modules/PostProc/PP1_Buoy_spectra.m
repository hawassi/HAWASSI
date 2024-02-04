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

idCol=1;
CR={'r','--g','-.b',':m','or','sg','db','+m'};
leg=[];
if EspId~=2
    SpectSaved=zeros(length(timeNow(1:floor(end/2))),length(buoyind));%% size ww==length(timeNow(2:floor(end/2));
    sppeakindexPP=0;
    nFreqP=str2num(dynmodel(3));
    BandFilt=GUIpp.Filter;
    if BandFilt==1
        LF=GUIpp.LFHF(1);      %Low Frequency
        HF=GUIpp.LFHF(2);     %High Frequency
    end
    
    for i=1:length(buoyind)
          
        if BandFilt==1
        simul=bandpass(timeNow,eta(indti:indtf,buoyind(i)),LF,HF);
        else
        simul=eta(indti:indtf,buoyind(i));
        end
        
        
        [ww,spSIMUL] = spect(timeNow,simul,Spsmooth,EspId);
        Nw=length(ww);indEn=floor(Nw); %6 minimal cutfracwn
        if strcmp(wall.presence,'Yes')
         
            if wall.type==2
                Omk=wall.Om;refl_coef=wall.refl_Coef;
                indOmend=floor(length(Omk)/2);
%                 par1=evalin('base','parR1');
%                 refl_coef1=par1.wall.refl_Coef;
%                 
                refl_coefN=interp1(Omk(1:indOmend),refl_coef(1:indOmend),ww(1:indEn),'spline');
                refl_coefN(isnan(refl_coefN))=refl_coef(indOmend);
%                 refl_coefN1=interp1(Omk(1:indOmend),refl_coef1(1:indOmend),ww(1:indEn),'spline');
%                 refl_coefN1(isnan(refl_coefN1))=refl_coef1(indOmend);
                
                if indOmend<indEn
                    refl_coefN(indOmend+1:indEn)=refl_coef(indOmend);
                end

                [AX,H1,H2]=plotyy(ww(1:dplot:indEn), spSIMUL(1:dplot:indEn),...
                    ww(1:dplot:indEn),refl_coefN(1:dplot:indEn));
%                  hold(AX(2))
%                  plot(AX(2),ww(1:dplot:indEn),refl_coefN1(1:dplot:indEn))
                
            else
                if idCol==9, idCol=1; end;
                plot(ww(1:dplot:indEn), spSIMUL(1:dplot:indEn),CR{idCol});
              
                if EspId==0,ylabel('Amplitude Sp.');
                else ylabel('Variance dens. Sp.');end;
                idCol=idCol+1;
                leg{i}=['Sp. x=',num2str(Xbuoy(i)),'[m]'];
                
            end
        else
            if idCol==9, idCol=1; end;
            plot(ww(1:dplot:indEn), spSIMUL(1:dplot:indEn),CR{idCol});
            if EspId==0,ylabel('Amplitude Sp.');
            else ylabel('Variance dens. Sp.');end;
            idCol=idCol+1;
            leg{i}=['Sp. x=',num2str(Xbuoy(i)),'[m]'];
        end
        hold on;
        SpectSaved(:,i)=spSIMUL;
        SpectTemp=abs(fft(simul));SpectTemp(1)=0;
        [~,IndexMaxSp]=max(SpectTemp);
        
        if IndexMaxSp>sppeakindexPP
           sppeakindexPP= IndexMaxSp;
        end
    end
    [MMsp,~]=max(max(SpectSaved));
    
    
    PPdata.PP1_Spectrum=SpectSaved;
    PPdata.PP1_frequency_rad=ww;
    if EspId ==1
        PPdata.PP1_Spectrum_type='Variance dens. Sp.';
    else
        PPdata.PP1_Spectrum_type='Amplitude';
    end
    
    if idCol~=2
       legend(leg{1:end});   
    end
    
    xlabel('\omega [rad/s]');
    
end

%nupeakPP=ww(sppeakindexPP);
nupeakPP=max(influx.nu_p,ww(sppeakindexPP));

if EspId~=2
    if strcmp(wall.presence,'Yes')
        if wall.type==2
            if nupeakPP~=0
                set(AX(1),'XLim',[0 min(nupeakPP*7,ww(indEn))])
                set(AX(2),'XLim',[0 min(nupeakPP*7,ww(indEn))])
            else
                set(AX(1),'XLim',[0 ww(indEn)/6])
                set(AX(2),'XLim',[0 ww(indEn)/6])
            end
            set(AX(2),'YLim',[0 1])
            set(AX(2),'Ytick',[0 0.5 1])
            if EspId==0
                ylabel(AX(1), 'Amplitude Sp.');
            else
                ylabel(AX(1), 'Variance dens. Sp.');
            end
            
            ylabel(AX(2), 'Reflection Coef.');
            
        else
            if nupeakPP~=0
                xlim([0 min(nupeakPP*7,ww(indEn))]);
            else
                xlim([0 ww(indEn)/6]);
            end
        end
    else
        if nupeakPP~=0
            xlim([0 min(nupeakPP*7,ww(indEn))]);
        else
            xlim([0 ww(indEn)/6]);
        end
    end
else
    if nupeakPP~=0
        xlim([0 min(nupeakPP*7,ww(indEn))]);
    else
        xlim([0 ww(indEn)/6]);
    end
end
% ylim([0 MMsp+MMsp/4]);
% MMspp=roundn(MMsp/4,-2);
% 
% if MMspp~=0
%     set(gca,'YTick',[0 2*MMspp 4*MMspp]);
% end
Nx=length(Xbuoy);
if Nx~=0
    if Nx~=1
        for i=1:Nx-1
            Xbuoy_str{i}=[num2str(Xbuoy(i)),';',];
        end
    end
    Xbuoy_str{Nx}=num2str(Xbuoy(Nx));
end;
title([dynmodel,',  Spectrum @ position: ',[Xbuoy_str{1:Nx}], ' [m]']);

% if GUIpp.PPsetting_ylim==1
%         ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
% end
if GUIpp.PPsetting_xlim==1
        xlim([GUIpp.PPsetting_xlimVal(1) GUIpp.PPsetting_xlimVal(2)])
end
if GUIpp.PPsetting_ylim==1
        ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
end
    plot_properties;