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

if GUIpp.PPsetting_tcoarse==1
    dplot= GUIpp.PPsetting_tcoarseVal;
else
    dplot=1;
end

Xbuoy   =GUIpp.PP1_Buoy_X';
buoyind = closest(x,Xbuoy);
etabuoy = eta(:,buoyind);

PPdata.PP1.Buoy_x=Xbuoy;
PPdata.PP1.Buoy_eta =etabuoy ;
PPdata.PP1.Buoy_time=time(1:dplot:end);

idCol=1;
CR={'r','--g','-.b',':m','or','sg','db','+m'};
leg=[];
if GUIpp.PP1_Buoy_Signal==1 && GUIpp.PP1_Buoy_Spectrum==0
    if GUIpp.PP1_Prof==1 && GUIpp.PP1_Buoy==1
        subplot(2,1,2)
    end
    
    if  strcmp(dynmodel,'HS1FWall')||strcmp(dynmodel,'HS2FWall')||strcmp(dynmodel,'HS3FWall')||strcmp(dynmodel,'HS4FWall')
        if isempty(handles.input)  || Load_data_Id==1;
            par=GUIpp.MySim.par;
            influx_dat=GUIpp.MySim.influx;
        else
            par=evalin('base','par');
            influx_dat=evalin('base','influx');
        end
       
        if isfield(par,'wall')
            reflPercent=par.wall.refl_Coef*100;
            Xwall=par.wall.position;
        else
            reflPercent=par.reflPercent;
            Xwall=par.wallposition;
        end
        
%         kp=influx_dat.k_p;
%         ampli=influx_dat.Hs./2;
%         TT=influx_dat.Tp;
%         Xinflux=influx_dat.position;
%         Xwall=x(closest(x,Xwall));
%         Xinflux=x(closest(x,Xinflux)-1);
%         lambda_nu=2*pi/kp;
%         hold on;
%         x_inp=Xbuoy;dx=x(2)-x(1);tsnap=PPdata.PP1.Buoy_time;
%         Eta_Sol1=ampli.*cos(kp.*(x_inp-Xinflux)-tsnap.*2*pi./TT);%+ampli.*cos(2*Xwall+Xinflux-xinp-t_inp.*2*pi/TT);
%         Eta_Sol2=(reflPercent/100).*ampli.*cos(-kp.*(x_inp-(Xinflux+2*(Xwall-Xinflux)))-tsnap.*2*pi./TT);%+ampli.*cos(2*Xwall+Xinflux-xinp-t_inp.*2*pi/TT);
%         Eta_Sol=Eta_Sol1+Eta_Sol2;
%         plot(tsnap,Eta_Sol,'--b')
        
    end
      BandFilt=GUIpp.Filter;
    if BandFilt==1
        LF=GUIpp.LFHF(1);      %Low Frequency
        HF=GUIpp.LFHF(2);     %High Frequency
    end
    
    for I=1:length(Xbuoy)
        if mod(idCol,9)==0, idCol=1; end
        if BandFilt==1
        etaFilt=bandpass(time,etabuoy(:,I),LF,HF);
        etaPlot=etaFilt(1:dplot:end);
        else
         etaPlot=etabuoy(1:dplot:end,I);   
        end
        plot(time(1:dplot:end), etaPlot,CR{idCol});
        hold on;
        idCol=idCol+1;
        leg{I}=['\eta x=',num2str(Xbuoy(I)),'[m]'];
    end
    
    xlabel('time [s]') ;   ylabel('\eta [m]');
    if idCol~=2
        legend(leg{1:end},'Location','northwest');
    end
    if GUIpp.PPsetting_xlim==1
        xlim([GUIpp.PPsetting_xlimVal(1) GUIpp.PPsetting_xlimVal(2)])
    else
        xlim([tinterv(1) tinterv(2)])
    end
    
    if GUIpp.PPsetting_ylim==1
        ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
    end
            
    plot_properties;
    if strcmp(wall.presence,'Yes')
    else
        Nx=length(Xbuoy);
        if Nx~=0
            if Nx~=1
                for i=1:Nx-1
                    Xbuoy_str{i}=[num2str(Xbuoy(i)),';',];
                end
            end
            Xbuoy_str{Nx}=num2str(Xbuoy(Nx));
        end;
        title([dynmodel,', Time signal @ position: ',[Xbuoy_str{1:Nx}], ' [m]']);
    end
    
elseif GUIpp.PP1_Buoy_Signal==0 && GUIpp.PP1_Buoy_Spectrum==1
    if GUIpp.PP1_Prof==1 && GUIpp.PP1_Buoy==1
        subplot(2,1,2)
    end
    
    
    if GUIpp.PP1_Buoy_Spectrum_Energy==1 &&GUIpp.PP1_Buoy_Spectrum_Ampl==0
        EspId=1;
    elseif GUIpp.PP1_Buoy_Spectrum_Ampl==1 && GUIpp.PP1_Buoy_Spectrum_Energy==0
        EspId=0;
    else
        EspId=2;
    end
    
    Spsmooth=GUIpp.PP1_Buoy_Spectrum_Smooth;
    
    indti=closest(time,tinterv(1));indtf=closest(time,tinterv(2));
    timeNow=time(indti:indtf);
    
    PP1_Buoy_spectra();
    
    
elseif GUIpp.PP1_Buoy_Signal==1 && GUIpp.PP1_Buoy_Spectrum==1
    if GUIpp.PP1_Prof==0 && GUIpp.PP1_Buoy==1
        subplot(2,1,1)
    elseif GUIpp.PP1_Prof==1 && GUIpp.PP1_Buoy==1
        subplot(3,1,2)
    end
    
      BandFilt=GUIpp.Filter;
    if BandFilt==1
        LF=GUIpp.LFHF(1);      %Low Frequency
        HF=GUIpp.LFHF(2);     %High Frequency
    end
    
    for I=1:length(Xbuoy)
        if idCol==9, idCol=1; end;
        if BandFilt==1
        etaFilt=bandpass(time,etabuoy(:,I),LF,HF);
        etaPlot=etaFilt(1:dplot:end);
        else
         etaPlot=etabuoy(1:dplot:end,I);   
        end
        
        plot(time(1:dplot:end), etaPlot,CR{idCol});
        idCol=idCol+1;
        hold on;
        leg{I}=['\eta x=',num2str(Xbuoy(I)),'[m]'];
    end
    xlabel('time [s]') ;   ylabel('\eta [m]');
    if idCol~=2
        legend(leg{1:end},'Location','northwest');
    end
    
    if GUIpp.PPsetting_xlim==1
        xlim([GUIpp.PPsetting_xlimVal(1) GUIpp.PPsetting_xlimVal(2)])
    else
        xlim([tinterv(1) tinterv(2)])
    end
    if GUIpp.PPsetting_ylim==1
        ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
    end
    plot_properties;
    
    Nx=length(Xbuoy);
    if Nx~=0
        if Nx~=1
            for i=1:Nx-1
                Xbuoy_str{i}=[num2str(Xbuoy(i)),';',];
            end
        end
        Xbuoy_str{Nx}=num2str(Xbuoy(Nx));
    end
    title([dynmodel,', Time signal @ position: ',[Xbuoy_str{1:Nx}], ' [m]']);
    
    
    if GUIpp.PP1_Prof==0 && GUIpp.PP1_Buoy==1
        subplot(2,1,2)
    elseif GUIpp.PP1_Prof==1 && GUIpp.PP1_Buoy==1
        subplot(3,1,3)
    end
    
    if GUIpp.PP1_Buoy_Spectrum_Energy==1 &&GUIpp.PP1_Buoy_Spectrum_Ampl==0
        EspId=1;
    elseif GUIpp.PP1_Buoy_Spectrum_Ampl==1 && GUIpp.PP1_Buoy_Spectrum_Energy==0
        EspId=0;
    else
        EspId=2;
    end
    
    Spsmooth=GUIpp.PP1_Buoy_Spectrum_Smooth;
    
    indti=closest(time,tinterv(1));indtf=closest(time,tinterv(2));
    timeNow=time(indti:indtf);
    
    PP1_Buoy_spectra();
end

if GUIpp.PP1_Quant==1
    PP1_Buoy_quant;
end

if GUIpp.PP1_Buoy_exceedance==1
   PP1_Buoy_exceedance; 
end