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

%%%%%% Post processing %%%%%%%%%%%%%%
global PPdata
if GUIpp.PP1==1,PP=1;
elseif GUIpp.PP2==1,PP=2;
elseif GUIpp.PP3==1,PP=3;
end


if isempty(handles.input)  || Load_data_Id==1
    if isfield(GUIpp.MySim,'model') %just indicates from the new version of saved data
        time=GUIpp.MySim.output.time;
        eta=GUIpp.MySim.output.eta;
        u=GUIpp.MySim.output.u;
        x=GUIpp.MySim.output.x;
        Xinflux=GUIpp.MySim.par.Xinflux;
        bath   =GUIpp.MySim.bath;
        influx  =GUIpp.MySim.influx;
        Oprt   =GUIpp.MySim.Oprt;
        bathy=GUIpp.MySim.par.bathy;
        minbath=-min(bathy);
        dynmodel=GUIpp.MySim.model.dyn;
        if isfield(GUIpp.MySim.par,'interp')
        OprtInterp_par=GUIpp.MySim.par.interp;
        end
        
        wall=GUIpp.MySim.par.wall;
        lambda_p=GUIpp.MySim.influx.lambda_p; % it is used  in the plot of Kinematic
        Nonlin_Adj=GUIpp.MySim.bath.influx_AdjZone;
        
        if isfield(GUIpp.MySim,'shippar')
            shippar=GUIpp.MySim.shippar;
            if shippar.check==1
                shipRB= GUIpp.MySim.output.RB;
                shipsavevar= GUIpp.MySim.output.shipsavevar;
            end
        else shippar.check=0;
        end
        
    else
        time=GUIpp.MySim.time;
        eta=GUIpp.MySim.eta;
        x=GUIpp.MySim.x;
        Xinflux=GUIpp.MySim.Xinflux;
        bathy=GUIpp.MySim.bathy;
        minbath=-min(bathy);
        dynmodel=GUIpp.MySim.dynmodel;
    end
    
    
    
    
    % if isdeployed
    
    
    if isfield(GUIpp.MySim,'dataBreak_nodes')
        dataBreak_nodes=GUIpp.MySim.dataBreak_nodes;
    else
        dataBreak_nodes=[];
    end
    
    if isfield(GUIpp.MySim,'output')
        dataBreak_nodes=GUIpp.MySim.output.break_nodes;
    end
    
    if isfield(GUIpp.MySim.output,'break_crest')
        break_crest=GUIpp.MySim.output.break_crest;
    else
        break_crest=[];
    end
    
    if PP==3
        if isfield(GUIpp.MySim,'model') %just indicates from the new version of saved data
            if isfield(GUIpp.MySim.influx,'tcoarse')
                tcoarse=GUIpp.MySim.influx.tcoarse;
            else
                tcoarse=1;
            end
        else
            if isfield(GUIpp.MySim,'tcoarse')
                tcoarse=GUIpp.MySim.tcoarse;
            else
                tcoarse=1;
            end
        end
        
        if isfield(GUIpp.MySim.par,'bf0')
            Fbdy=GUIpp.MySim.par.bf0;
        else
            Fbdy=0;
        end
    end
    
elseif ~isempty(handles.input)
    
    time=handles.input{1}.time;
    eta=handles.input{1}.eta;
    u=handles.input{1}.u;
    x=handles.input{1}.x;
    bathy=handles.input{1}.bathy;
    bath   =handles.input{1}.bath;
    Oprt   =handles.input{1}.Oprt;
    dynmodel=handles.input{1}.dynmodel;
    Xinflux=handles.input{1}.Xinflux;
    influx=handles.input{1}.influx;
    if isfield(handles.input{1},'OprtInterp_par')
        OprtInterp_par=handles.input{1}.OprtInterp_par;
    end
    
    minbath=-min(bathy);
    wall=handles.input{1}.Wall;
    dataBreak_nodes=handles.input{1}.dataBreak_nodes;
    break_crest=handles.input{1}.break_crest;
    lambda_p=handles.input{1}.lambda_p; % it is used  in the plot of Kinematic
    Nonlin_Adj=handles.input{1}.Nonlin_Adj;
    shippar=handles.input{1}.shippar;
    if shippar.check==1
        shipRB=handles.input{1}.shipRB;
        shipsavevar= handles.input{1}.shipsavevar;
    end
    
    if PP==3
        %         IndtI=handles.input{1}.IndtI;
        %         IndtF=handles.input{1}.IndtF;
        if isfield(handles.input{1},'tcoarse')
            tcoarse=handles.input{1}.tcoarse;
        else
            tcoarse=1;
        end
        Fbdy=handles.input{1}.Fbdy;
    end
    
end

if isempty(dataBreak_nodes)
    if GUIpp.PP1_breakevents==1 &&GUIpp.PP1==1
        statusbarObj.setText('No breaking data available!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return;
    end
end

sf_savename         = GUIpp.workdir;

%%%%%%%%%%To keep adding data during PP process%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%it will change if new file is loaded%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp(GUIpp.workdir,PPdata.workdir)
    PPdata=[];
    PPdata.workdir=GUIpp.workdir;
    PPdata.PP1.name='Profile and time signals';
    PPdata.PP2.name='Animation: No data saved';
    PPdata.PP3.name='Validation';
end



if PP==1
    if GUIpp.PPsetting_MTA_time==1
        indt1=closest(time,GUIpp.PPsetting_MTA_timeVal(1));
        indt2=closest(time,GUIpp.PPsetting_MTA_timeVal(2));
        MTC     = max(eta(indt1:indt2,:));  MTT     = min(eta(indt1:indt2,:));
    else
        MTC     = max(eta);  MTT     = min(eta);
    end
    % Hs      = 4*sqrt(var(eta,1));
    
    if  GUIpp.PP1_Prof==1||GUIpp.PP1_Buoy==1
        if   GUIpp.PP1_Buoy==1
            if GUIpp.PP1_Buoy_Signal==1 ||GUIpp.PP1_Buoy_Spectrum==1
                hf1=figure('Name','PostProc', 'Position',[250,10,800,500]);
            end
        else
            hf1=figure('Name','PostProc', 'Position',[250,10,800,500]);
        end
    end
    IdF=0;
    if  GUIpp.PP1_Prof==1 && GUIpp.PP1_Buoy==0
        PP1_Profile();IdF=1;
    elseif  GUIpp.PP1_Prof==0 && GUIpp.PP1_Buoy==1
        PP1_Buoy();IdF=1;
    elseif   GUIpp.PP1_Prof==1 && GUIpp.PP1_Buoy==1
        if (GUIpp.PP1_Buoy_Signal==1 &&GUIpp.PP1_Buoy_Spectrum==0) ...
                || (GUIpp.PP1_Buoy_Signal==0 &&GUIpp.PP1_Buoy_Spectrum==1)
            subplot(2,1,1)
        elseif (GUIpp.PP1_Buoy_Signal==1 &&GUIpp.PP1_Buoy_Spectrum==1)
            subplot(3,1,1)
        end
        PP1_Profile();
        PP1_Buoy();
        IdF=1;
    end
    
    if GUIpp.PP1_Spectra==1
       PP1_Spectra; 
    end
    
    if GUIpp.PP1_Energy==1
        PP1_Energy;
    end
    
    if GUIpp.PP1_Energy_in_x==1
       PP1_Energy_in_x; 
    end
    
    if GUIpp.PP1_kinematic==1
        PP1_Kinematic;
    end
    
    if GUIpp.PP1_breakevents==1
        if isempty(break_crest) && ~isempty(dataBreak_nodes)
            PP1_BreakEvents;
        else
            hf3=figure('Name','PostProc', 'Position',[250,10,800,500]);
            tbreak=break_crest(:,1);
            Nodes=break_crest(:,2:end);
            Nodes(Nodes==0)=1;
            xb=x(Nodes);
            xb(Nodes==1)=NaN;
            plot(tbreak,xb,'ob');
            
            ylabel('Breaking position [m]');
            xlabel('time [s]');
            PPdata.PP1.breaking.position=xb;
            PPdata.PP1.breaking.time=tbreak;
            plot_properties;
        end
    end
    
    if GUIpp.PP1_MTA_MWL==1
        PP1_MTA_MWL;
    end
    
    if GUIpp.PP1_HsAsSkK==1
        PP1_Hs_Sk_As_K;
    end
    
    if GUIpp.PP1_Ship_XZThe==1
        if shippar.check==1
        PP1_Ship_XZThe;
        else
        statusbarObj.setText('Ship data is not available!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);  
        end
    end
    
    if GUIpp.PP1_Ship_Forces==1 && IdF==0
        if shippar.check==1
            PP1_Ship_Forces;
        else
            statusbarObj.setText('Ship data is not available!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
        end
    end
    
    
    if GUIpp.save_figure==1
        if  GUIpp.PP1_Prof==1||GUIpp.PP1_Buoy==1
            if strcmp(GUIpp.save_figure_type,'.eps')
                saveas(hf1,[sf_savename,'PP1_Profile_Signal',GUIpp.save_figure_type],'epsc')
            else
                saveas(hf1,[sf_savename,'PP1_Profile_Signal',GUIpp.save_figure_type])
            end
        end
        if GUIpp.PP1_Energy==1
            if strcmp(GUIpp.save_figure_type,'.eps')
                saveas(hf2,[sf_savename,'PP1_Energy',GUIpp.save_figure_type],'epsc')
            else
                saveas(hf2,[sf_savename,'PP1_Energy',GUIpp.save_figure_type])
            end
        end
        if GUIpp.PP1_breakevents==1
            if strcmp(GUIpp.save_figure_type,'.eps')
                saveas(hf3,[sf_savename,'PP1_breakevents',GUIpp.save_figure_type],'epsc')
            else
                saveas(hf3,[sf_savename,'PP1_breakevents',GUIpp.save_figure_type])
            end
        end
        if GUIpp.PP1_kinematic==1
            if strcmp(GUIpp.save_figure_type,'.eps')
                saveas(hf4,[sf_savename,'PP1_Kinematic',GUIpp.save_figure_type],'epsc')
            else
                saveas(hf4,[sf_savename,'PP1_Kinematic',GUIpp.save_figure_type])
            end
        end
        
        if GUIpp.PP1_Ship_XZThe==1
            if strcmp(GUIpp.save_figure_type,'.eps')
                saveas(hf5,[sf_savename,'PP1_Ship_XZTheta',GUIpp.save_figure_type],'epsc')
            else
                saveas(hf5,[sf_savename,'PP1_Ship_XZTheta',GUIpp.save_figure_type])
            end
        end
        
        if GUIpp.PP1_Ship_Forces==1 && IdF==0
            if strcmp(GUIpp.save_figure_type,'.eps')
                saveas(hf6,[sf_savename,'PP1_Ship_Forces',GUIpp.save_figure_type],'epsc')
            else
                saveas(hf6,[sf_savename,'PP1_Ship_Forces',GUIpp.save_figure_type])
            end
        end
        
    end
    
    %
    %     figure
    %     plot(x,Eta_Sol1,'b',x,Eta_Sol2,'--r',x,Eta_Sol,'-c')
    
    
elseif PP==2
    BandFilt=GUIpp.Filter;
    if BandFilt==1
        statusbarObj.setText('Filtering data...');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.black);
            
        LF=GUIpp.LFHF(1);      %Low Frequency
        HF=GUIpp.LFHF(2);      %High Frequency
        etaFilt=zeros(size(eta));
        for ii=1:length(x)
            etaFilt(:,ii)=bandpass(time,eta(:,ii),LF,HF);
        end
        eta=etaFilt;
        statusbarObj.setText('');
        statusbarTxt = statusbarObj.getComponent(0);
    end
   
            
            
    
    
    if GUIpp.PPsetting_MTA_time==1
        indt1=closest(time,GUIpp.PPsetting_MTA_timeVal(1));
        indt2=closest(time,GUIpp.PPsetting_MTA_timeVal(2));
        MTC     = max(eta(indt1:indt2,:));  MTT     = min(eta(indt1:indt2,:));
    else
        MTC     = max(eta);  MTT     = min(eta);
    end
    if shippar.check==1
        for jj=1:shippar.Nship
            IndSL=closest(x,shippar.form.xShip(jj,1));
            IndSR=closest(x,shippar.form.xShip(jj,3));
            
            MTT(IndSL+1:IndSR-1)=NaN;
            MTC(IndSL+1:IndSR-1)=NaN; 
%             figure;
%             plot(x,MTT,'r',x,MTC,'b')
%       
      
            draft=shippar.form.Sdraft(jj);
            slength=shippar.form.Slength(jj);
            
            z0SS=linspace(-draft,1.2*max(MTC),20);
            x0SS=linspace(-slength/2,slength/2,20);
            onez=ones(size(z0SS));
            onex=ones(size(x0SS));
            
            S(jj).zSS=[z0SS       -draft.*onex    z0SS            1.2*max(MTC).*onex];
            S(jj).xSS=[x0SS(1)*onez x0SS          x0SS(end).*onez fliplr(x0SS)];
            S(jj).xS0=shippar.form.xShip(jj,2);
        end 
          
          
%             figure;
%             fill(xSSn0,zSSn0,'r')
    end
        
    
    scale=GUIpp.PP2_Bathy_scale;
    GIF_delaytime=GUIpp.PP2_GIF_setting(1);
    GIF_loopcount=GUIpp.PP2_GIF_setting(2);
    
    indti=closest(time,GUIpp.PP2_T(1));
    indtf=closest(time,GUIpp.PP2_T(2));
    
    dt   =time(2)-time(1);
    
    
    if GUIpp.PPsetting_tcoarse==1
        step= GUIpp.PPsetting_tcoarseVal;
    else
        step=1;                              %round(GUIpp.PP4.T(2)/(time(2)-time(1)));
    end
    
    xinit=GUIpp.PP2_X(1);xend=GUIpp.PP2_X(2);
    
    
    if GUIpp.PPsetting_xcoarse==1
        stepdx= GUIpp.PPsetting_xcoarseVal;
    else
        stepdx=1;                              %round(GUIpp.PP4.T(2)/(time(2)-time(1)));
    end
    
    
    indxI=1;indxE=length(x);
    filename = [sf_savename,'PP2_animation.gif'];
    assignin('base','sf_savename',sf_savename);
    
    if ~isfolder(sf_savename)
        mkdir(sf_savename);
    end
    hf2=figure;
   set(hf2, 'units','normalized','Position', [0.2, 0.2, 0.6,0.5]);
    %set(hf2, 'units','normalized','Position', [0.05, 0.1, 1,0.8]);
    
%     ax = gca;
%     outerpos = ax.OuterPosition;
%     ti = ax.TightInset;
%     left = outerpos(1) + ti(1);
%     bottom = outerpos(2) + ti(2);
%     ax_width = outerpos(3) - ti(1) - ti(3);
%     ax_height = outerpos(4) - ti(2) - ti(4);

         
    numframes=length(time(indti:step:indtf));
    MM = roundn(max(max(eta)),-1);
    
    
    bathy_plot=bathy*scale;
    
    BreakModID=0;
    if ~isempty(dataBreak_nodes), BreakModID=1;end;
    if  strcmp(wall.presence,'Yes')
        ampli=max(max(eta));
        Ywall=(-2*ampli:ampli/5:2*ampli);
        Xwall=wall.position*ones(size(Ywall));
        
    end
    
    
    
    for i = 1:numframes
        figure(hf2)
        indt=indti+i*step-1;
        %             H=eta(indt,:)+hb;
        %             HeavH=Heaviside(H-H_min);
        %             ind0=find(HeavH==0,1,'first'); %shorepoint index
        %             ind1=ind0-1;
        if indt>indtf, break; end
        
        if shippar.check==1
            Nship=shippar.Nship;
             sX=real(shipRB(indt,1:Nship));
             for jjNs=1:Nship
                 IndSL=closest(x,shippar.form.xShip(jjNs,1)+sX(jjNs));
                 IndSR=closest(x,shippar.form.xShip(jjNs,3)+sX(jjNs));
                 eta(indt,IndSL+1:IndSR-1)=nan;
             end
        end
        
        if GUIpp.PP2_MTA==0 && GUIpp.PP2_Bathy==1
            
            if  strcmp(wall.presence,'Yes')
                if BreakModID==1&&time(indt)>=dataBreak_nodes(1,1)&&time(indt)<=dataBreak_nodes(end,1)
                    indtBreak=closest(dataBreak_nodes(:,1),time(indt));
                    indendb=find(dataBreak_nodes(indtBreak,:)==0,1,'first')-1;
                    indxBreak=dataBreak_nodes(indtBreak,2:indendb);
                    
                    if all(abs(time(indt)-dataBreak_nodes(:,1))>dt)
                        plot (x(indxI:stepdx:indxE),eta(indt,indxI:stepdx:indxE),'r',x(indxI:stepdx:indxE),bathy_plot(indxI:stepdx:indxE),'--k',Xwall,Ywall,'oy');
                    else
                        plot (x(indxI:stepdx:indxE),eta(indt,indxI:stepdx:indxE),'r',x([indxBreak]),eta(indt,[indxBreak]),'ob', ...
                            Xinflux,eta(indti+i*step-1,closest(x,Xinflux)),'c*',x(indxI:indxE),bathy_plot(indxI:indxE),'k--',Xwall,Ywall,'oy');
                    end
                else
                    plot (x(indxI:stepdx:indxE),eta(indt,indxI:stepdx:indxE),'r',x(indxI:stepdx:indxE),bathy_plot(indxI:stepdx:indxE),'--k',Xwall,Ywall,'oy');
                end
            else
                if BreakModID==1&&time(indt)>=dataBreak_nodes(1,1)&&time(indt)<=dataBreak_nodes(end,1)
                    indtBreak=closest(dataBreak_nodes(:,1),time(indt));
                    indendb=find(dataBreak_nodes(indtBreak,:)==0,1,'first')-1;
                    indxBreak=dataBreak_nodes(indtBreak,2:indendb);
                    
                    if all(abs(time(indt)-dataBreak_nodes(:,1))>dt)
                        plot (x(indxI:stepdx:indxE),eta(indt,indxI:stepdx:indxE),'r',x(indxI:stepdx:indxE),bathy_plot(indxI:stepdx:indxE),'--k');
                    else
                        plot (x(indxI:stepdx:indxE),eta(indt,indxI:stepdx:indxE),'r',x([indxBreak]),eta(indt,[indxBreak]),'ob', ...
                            Xinflux,eta(indti+i*step-1,closest(x,Xinflux)),'c*',x(indxI:indxE),bathy_plot(indxI:indxE),'k--');
                    end
                else
                    plot (x(indxI:stepdx:indxE),eta(indt,indxI:stepdx:indxE),'r',x(indxI:stepdx:indxE),bathy_plot(indxI:stepdx:indxE),'--k');
                end 
            end
            
            if any(bathy>0)
                ylim([min(min(MTT),min(bathy_plot)) max(MTC(bathy<=0))]);
            else
                ylim([min(min(MTT),min(bathy_plot)) max(MTC)]);
            end
        elseif GUIpp.PP2_MTA==1 && GUIpp.PP2_Bathy==0
            if  strcmp(wall.presence,'Yes')
                if BreakModID==1&&time(indt)>=dataBreak_nodes(1,1)&&time(indt)<=dataBreak_nodes(end,1)
                    indtBreak=closest(dataBreak_nodes(:,1),time(indt));
                    indendb=find(dataBreak_nodes(indtBreak,:)==0,1,'first')-1;
                    indxBreak=dataBreak_nodes(indtBreak,2:indendb);
                    
                    if all(abs(time(indt)-dataBreak_nodes(:,1))>dt)
                        plot (x(indxI:stepdx:indxE),eta(indt,indxI:stepdx:indxE),'r',...
                            x(indxI:stepdx:indxE),MTC(indxI:stepdx:indxE),'-.k',...
                            x(indxI:stepdx:indxE),MTT(indxI:stepdx:indxE),':c',...
                            Xwall,Ywall,'oy');
                    else
                        plot (x(indxI:stepdx:indxE),eta(indt,indxI:stepdx:indxE),'r',...
                            x([indxBreak]),eta(indt,[indxBreak]),'ob', ...
                            Xinflux,eta(indti+i*step-1,closest(x,Xinflux)),'c*',...
                            x(indxI:stepdx:indxE),MTC(indxI:stepdx:indxE),'-.k',...
                            x(indxI:stepdx:indxE),MTT(indxI:stepdx:indxE),':c',...
                            Xwall,Ywall,'oy');
                    end
                else
                    plot (x(indxI:stepdx:indxE),eta(indt,indxI:stepdx:indxE),'r',...
                        x(indxI:stepdx:indxE),MTC(indxI:stepdx:indxE),'-.k',...
                        x(indxI:stepdx:indxE),MTT(indxI:stepdx:indxE),':c',...
                        Xwall,Ywall,'oy');
                end
            else
                if BreakModID==1&&time(indt)>=dataBreak_nodes(1,1)&&time(indt)<=dataBreak_nodes(end,1)
                    indtBreak=closest(dataBreak_nodes(:,1),time(indt));
                    indendb=find(dataBreak_nodes(indtBreak,:)==0,1,'first')-1;
                    indxBreak=dataBreak_nodes(indtBreak,2:indendb);
                    
                    if all(abs(time(indt)-dataBreak_nodes(:,1))>dt)
                        plot (x(indxI:stepdx:indxE),eta(indt,indxI:stepdx:indxE),'r',x(indxI:stepdx:indxE),MTC(indxI:stepdx:indxE),'-.k',x(indxI:stepdx:indxE),MTT(indxI:stepdx:indxE),':c');
                    else
                        plot (x(indxI:stepdx:indxE),eta(indt,indxI:stepdx:indxE),'r',x([indxBreak]),eta(indt,[indxBreak]),'ob', ...
                            Xinflux,eta(indti+i*step-1,closest(x,Xinflux)),'c*',...
                            x(indxI:stepdx:indxE),MTC(indxI:stepdx:indxE),'-.k',x(indxI:stepdx:indxE),MTT(indxI:stepdx:indxE),':c');
                    end
                else
                    plot (x(indxI:stepdx:indxE),eta(indt,indxI:stepdx:indxE),'r',x(indxI:stepdx:indxE),MTC(indxI:stepdx:indxE),'-.k',x(indxI:stepdx:indxE),MTT(indxI:stepdx:indxE),':c');
                    %                      indsl=closest(x,shippar.form.xShip(1));indsr=closest(x,shippar.form.xShip(3));
                    %                      etaS=eta(indt,:);MTCS=MTC;MTTS=MTT;
                    %                      etaS(indsl+1:indsr-1)=NaN;MTCS(indsl+1:indsr-1)=MTC(indsl+1:indsr-1)+12.5;MTTS(indsl+1:indsr-1)=MTT(indsl+1:indsr-1)+12.5;
                    %                      etaSS=eta(indt,indsl+1:indsr-1);
                    %                      etaSS=etaSS+12.5;
                    %
                    %                      plot (x(indxI:stepdx:indxE),etaS(indxI:stepdx:indxE),'r',x(indxI:stepdx:indxE),MTCS(indxI:stepdx:indxE),'-.k',x(indxI:stepdx:indxE),MTTS(indxI:stepdx:indxE),':c',...
                    %                            x(indsl+1:indsr-1),etaSS,'-ok');
                    
                    
                end
                
            end
            if any(bathy>0)
              ylim([min(MTT) max(MTC(bathy<=0))]);
            else
                if min(MTT)~=max(MTC)
                    ylim([min(MTT) max(MTC)]);
                end
            end
        elseif GUIpp.PP2_MTA==1 && GUIpp.PP2_Bathy==1
            if  strcmp(wall.presence,'Yes')
                if BreakModID==1&&time(indt)>=dataBreak_nodes(1,1)&&time(indt)<=dataBreak_nodes(end,1)
                    indtBreak=closest(dataBreak_nodes(:,1),time(indt));
                    indendb=find(dataBreak_nodes(indtBreak,:)==0,1,'first')-1;
                    indxBreak=dataBreak_nodes(indtBreak,2:indendb);
                    
                    if all(abs(time(indt)-dataBreak_nodes(:,1))>dt)
                        plot (x(indxI:stepdx:indxE),eta(indt,indxI:stepdx:indxE),'r',x(indxI:stepdx:indxE),MTC(indxI:stepdx:indxE),'-.k',x(indxI:stepdx:indxE),MTT(indxI:stepdx:indxE),':c',...
                            x(indxI:stepdx:indxE),bathy_plot(indxI:stepdx:indxE),'k--',Xwall,Ywall,'oy');
                    else
                        plot (x(indxI:stepdx:indxE),eta(indt,indxI:stepdx:indxE),'r',x([indxBreak]),eta(indt,[indxBreak]),'ob', ...
                            Xinflux,eta(indti+i*step-1,closest(x,Xinflux)),'c*',...
                            x(indxI:stepdx:indxE),MTC(indxI:stepdx:indxE),'-.k',x(indxI:stepdx:indxE),MTT(indxI:stepdx:indxE),':c',...
                            x(indxI:stepdx:indxE),bathy_plot(indxI:stepdx:indxE),'k--',Xwall,Ywall,'oy');
                    end
                else
                    plot (x(indxI:stepdx:indxE),eta(indt,indxI:stepdx:indxE),'r',x(indxI:stepdx:indxE),MTC(indxI:stepdx:indxE),'-.k',x(indxI:stepdx:indxE),MTT(indxI:stepdx:indxE),':c',...
                        x(indxI:stepdx:indxE),bathy_plot(indxI:stepdx:indxE),'k--',Xwall,Ywall,'oy');
                end
            else
                if BreakModID==1&&time(indt)>=dataBreak_nodes(1,1)&&time(indt)<=dataBreak_nodes(end,1)
                    indtBreak=closest(dataBreak_nodes(:,1),time(indt));
                    indendb=find(dataBreak_nodes(indtBreak,:)==0,1,'first')-1;
                    indxBreak=dataBreak_nodes(indtBreak,2:indendb);
                    
                    if all(abs(time(indt)-dataBreak_nodes(:,1))>dt)
                        plot (x(indxI:stepdx:indxE),eta(indt,indxI:stepdx:indxE),'r',x(indxI:stepdx:indxE),MTC(indxI:stepdx:indxE),'-.k',x(indxI:stepdx:indxE),MTT(indxI:stepdx:indxE),':c',...
                            x(indxI:stepdx:indxE),bathy_plot(indxI:stepdx:indxE),'k--');
                    else
                        plot (x(indxI:stepdx:indxE),eta(indt,indxI:stepdx:indxE),'r',x([indxBreak]),eta(indt,[indxBreak]),'ob', ...
                            Xinflux,eta(indti+i*step-1,closest(x,Xinflux)),'c*',...
                            x(indxI:stepdx:indxE),MTC(indxI:stepdx:indxE),'-.k',x(indxI:stepdx:indxE),MTT(indxI:stepdx:indxE),':c',...
                            x(indxI:stepdx:indxE),bathy_plot(indxI:stepdx:indxE),'k--');
                    end
                else
                    plot (x(indxI:stepdx:indxE),eta(indt,indxI:stepdx:indxE),'r',x(indxI:stepdx:indxE),MTC(indxI:stepdx:indxE),'-.k',x(indxI:stepdx:indxE),MTT(indxI:stepdx:indxE),':c',...
                        x(indxI:stepdx:indxE),bathy_plot(indxI:stepdx:indxE),'k--');
                end
                
            end
            
            if any(bathy>0)
                ylim([min(min(MTT),min(bathy_plot)) max(MTC(bathy<=0))]);
            else
                ylim([min(min(MTT),min(bathy_plot)) max(MTC)]);
            end
        elseif GUIpp.PP2_MTA==0 && GUIpp.PP2_Bathy==0
            if  strcmp(wall.presence,'Yes')
                if BreakModID==1&&time(indt)>=dataBreak_nodes(1,1)&&time(indt)<=dataBreak_nodes(end,1)
                    indtBreak=closest(dataBreak_nodes(:,1),time(indt));
                    indendb=find(dataBreak_nodes(indtBreak,:)==0,1,'first')-1;
                    indxBreak=dataBreak_nodes(indtBreak,2:indendb);
                    
                    if all(abs(time(indt)-dataBreak_nodes(:,1))>dt)
                        plot (x(indxI:stepdx:indxE),eta(indt,indxI:stepdx:indxE),'r',Xwall,Ywall,'oy');
                    else
                        plot (x(indxI:stepdx:indxE),eta(indt,indxI:stepdx:indxE),'r',x([indxBreak]),eta(indt,[indxBreak]),'ob', ...
                            Xinflux,eta(indti+i*step-1,closest(x,Xinflux)),'c*',Xwall,Ywall,'oy');
                    end
                else
                    plot (x(indxI:stepdx:indxE),eta(indt,indxI:stepdx:indxE),'r',Xwall,Ywall,'oy');
                end
            else
                if BreakModID==1&&time(indt)>=dataBreak_nodes(1,1)&&time(indt)<=dataBreak_nodes(end,1)
                    indtBreak=closest(dataBreak_nodes(:,1),time(indt));
                    indendb=find(dataBreak_nodes(indtBreak,:)==0,1,'first')-1;
                    indxBreak=dataBreak_nodes(indtBreak,2:indendb);
                    
                    if all(abs(time(indt)-dataBreak_nodes(:,1))>dt)
                        plot (x(indxI:stepdx:indxE),eta(indt,indxI:stepdx:indxE),'r');
                    else
                        plot (x(indxI:stepdx:indxE),eta(indt,indxI:stepdx:indxE),'r',x([indxBreak]),eta(indt,[indxBreak]),'ob', ...
                            Xinflux,eta(indti+i*step-1,closest(x,Xinflux)),'c*');
                    end
                else
                    plot (x(indxI:stepdx:indxE),eta(indt,indxI:stepdx:indxE),'r');
                end
                
            end
            
       
            
            
            if any(bathy>0)
                 ylim([min(MTT) max(MTC(bathy<=0))]);
            else
                ylim([min(MTT) max(MTC)]);
            end
        end
        
        if shippar.check==1
            sX=real(shipRB(indt,1:Nship));
            sZ=real(shipRB(indt,2*Nship+1:Nship*3));
            sTheta=real(shipRB(indt,4*Nship+1:5*Nship));
            for jj=1:Nship
                xSSn=S(jj).xS0+sX(jj)+(S(jj).xSS).*cos(sTheta(jj))-(S(jj).zSS).*sin(sTheta(jj));
                zSSn=sZ(jj)+(S(jj).xSS).*sin(sTheta(jj))+(S(jj).zSS).*cos(sTheta(jj));
                
                hold on
                fill(xSSn,zSSn,'b')
                hold off;
            end
      %    axis equal;
        end
                
        
        title(['t= ',num2str(time(indti+i*step-1)), ' [s]']);
        grid on; %   xlim([-10 10]);
        %         MM = ceil(max(max(eta)));
        if GUIpp.PPsetting_xlim==1
            xlim([GUIpp.PPsetting_xlimVal(1) GUIpp.PPsetting_xlimVal(2)])
        else
            xlim([xinit xend])
        end
        if GUIpp.PPsetting_ylim==1
            ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
        end
        
      
            
        xlabel('x[m]'); ylabel('\eta[m]')
        plot_properties;
        pause(GIF_delaytime)
        drawnow
        if ~ishandle(hf2)
            break;
        end
        
        if GUIpp.save_animation==1
            try
                frame =getframe(hf2);
            catch
                
            end
            
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            if i == 1;
                imwrite(imind,cm,filename,'gif','DelayTime',GIF_delaytime,'loopcount',GIF_loopcount);
            else
                imwrite(imind,cm,filename,'gif','DelayTime',GIF_delaytime, 'Writemode', 'append');
            end
        end
        
    end
  
    
elseif PP==3
    if GUIpp.PP3_AsSk==1
        JavProgressBar;
    end
    
    MeasPos=GUIpp.PP3_X;    %Measurement positions
    Nfig=length(MeasPos);        % Number of measurement
    
    Spsmooth=GUIpp.PP3_Spectra_smooth;   %Spectrum smoothing
    
    
    
    BandFilt=GUIpp.Filter;
    if BandFilt==1
        LF=GUIpp.LFHF(1);      %Low Frequency
        HF=GUIpp.LFHF(2);      %High Frequency
    end
    
    if GUIpp.PP3_Spectra_Energy==1 &&GUIpp.PP3_Spectra_Ampl==0
        ESpId=1;
    elseif GUIpp.PP3_Spectra_Ampl==1 && GUIpp.PP3_Spectra_Energy==0
        ESpId=0;
    elseif  GUIpp.PP3_Spectra_Ampl==1 && GUIpp.PP3_Spectra_Energy==1
        ESpId=2;
    else
        ESpId=1;
    end
    
    timemeas=GUIpp.PP3MeasDat(2:end,1);
    IndtI=closest(timemeas,time(1));
    IndtF=closest(timemeas,time(end));
    dtSim=time(2)-time(1);
    dtMeas=timemeas(2)-timemeas(1);
    
    MeasData=GUIpp.PP3MeasDat(1+IndtI:tcoarse:1+IndtF,:);   %Measurement data
    timemeasnow=timemeas(IndtI:tcoarse:IndtF);
    indti=closest(time,GUIpp.PP3_T(1));
    indtf=closest(time,GUIpp.PP3_T(2));
    
    if indtf>IndtF
        indtf=IndtF;
    end
    
    
    if GUIpp.PPsetting_tcoarse==1
        dplot= GUIpp.PPsetting_tcoarseVal;
    else
        dplot=1;                              %round(GUIpp.PP4.T(2)/(time(2)-time(1)));
    end
    
    if GUIpp.PPsetting_xcoarse==1
        stepdx= GUIpp.PPsetting_xcoarseVal;
    else
        stepdx=1;                              %round(GUIpp.PP4.T(2)/(time(2)-time(1)));
    end
    
    indtf=min(indtf,size(eta,1));
    indtf=min(indtf,size(MeasData,1));
    timeNow=time(indti:dplot:indtf);
    
    
    vari        = zeros(Nfig,3);
    Hs          = zeros(Nfig,3);
    Corr        = zeros(Nfig,1);CorrN=Corr;
    OptCorr     = zeros(Nfig,2);
    kurt        = zeros(Nfig,3);
    Asym        = zeros(Nfig,3);
    Skew        = zeros(Nfig,3);
    MWL         = zeros(Nfig,2);MTCbuoy=MWL;MTTbuoy=MWL;
    spSIMUL     = zeros(floor(length(timeNow)/2),Nfig);
    spMEAS      = zeros(floor(length(timeNow)/2),Nfig);
    
    
    if ESpId==2
        spSIMULEn     = zeros(floor(length(timeNow)/2),Nfig);
        spMEASEn      = zeros(floor(length(timeNow)/2),Nfig);
    end
    
    if GUIpp.PP3_Signals==1
        hf3a=figure;
    end
    FlagYlabel=0;
    
    simulSaved=zeros(length(indti:dplot:indtf),Nfig);
    measSaved=zeros(length(indti:dplot:indtf),Nfig);
    
    for i=1:Nfig
        indX=closest(x,MeasPos(i));
        IndXMeasdat=closest(GUIpp.PP3MeasDat(1,2:end),MeasPos(i))+1;
        if BandFilt==1
            etaFilt=bandpass(time,eta(:,indX),LF,HF);
            measFilt=bandpass(time,MeasData(1:length(time),IndXMeasdat),LF,HF);
            simul=etaFilt(indti:dplot:indtf);
            meas=measFilt(indti:dplot:indtf);
        else
            simul = eta(indti:dplot:indtf,indX);
            meas =  MeasData(indti:dplot:indtf,IndXMeasdat);
        end
        %         meas=meas-mean(meas);
        
        
        MM1 =max([max(simul) max(meas)]);
        MM =MM1+MM1/2;
        
        if ESpId~=2
            [w,spsimul,varsimul] = spect(timeNow,simul,Spsmooth,ESpId);
            [ww,spmeas,varmeas]  = spect(timeNow,meas,Spsmooth,ESpId);
        else
            [~,spsimul,~] = spect(timeNow,simul,Spsmooth,0);
            [~,spmeas,~]  = spect(timeNow,meas,Spsmooth,0);
            [w,spsimulEn,varsimul] = spect(timeNow,simul,Spsmooth,1);
            [ww,spmeasEn,varmeas]  = spect(timeNow,meas,Spsmooth,1);
            spSIMULEn(:,i)=spsimulEn'; spMEASEn(:,i)=spmeasEn';
        end
        spSIMUL(:,i)=spsimul'; spMEAS(:,i)=spmeas';
        
        
        vari(i,2)         = varmeas;
        vari(i,1)         = varsimul;
        vari(i,3)         = varsimul/varmeas;
        kurt(i,2)         = kurtosis(meas);
        kurt(i,1)         = kurtosis(simul);
        kurt(i,3)         = kurt(i,1)/kurt(i,2);
        Hs(i,2)   = 4*sqrt(varmeas);
        Hs(i,1)   = 4*sqrt(varsimul);
        Hs(i,3)   = Hs(i,1)/Hs(i,2);
        assignin('base','simul',simul);
        assignin('base','meas',meas);
        
       % Corr(i)   = (simul'*meas)/sqrt(varmeas*varsimul)/(length(meas));%correl(simul,meas);
        Corr(i)   = (simul'*meas)/norm(simul)/norm(meas);%correl(simul,meas);
       
        
        MWLsimul=mean(simul);
        MWLmeas =mean(meas);
        MWL(i,1)=MWLsimul;MWL(i,2)=MWLmeas;
        Hilb_meas=zeros(size(meas));%imag(hilbert(meas-MWLmeas));
        Hilb_sim=zeros(size(meas));%imag(hilbert(simul-MWLsimul));
        Asym(i,1)=mean(Hilb_sim.^3)/((mean((simul-MWLsimul).^2))^(3/2));
        Asym(i,2)=mean(Hilb_meas.^3)/((mean((meas-MWLmeas).^2))^(3/2));
        Asym(i,3)=Asym(i,1)/Asym(i,2);
        Skew(i,1)=mean((simul-MWLsimul).^3)/((mean((simul-MWLsimul).^2))^(3/2));
        Skew(i,2)=mean((meas-MWLmeas).^3)/((mean((meas-MWLmeas).^2))^(3/2));
        Skew(i,3)=Skew(i,1)/Skew(i,2);
        %Skew(i,1)=skewness(simul);
        %Skew(i,2)=skewness(meas);
        
        
        MTCbuoy(i,1)=max(simul);MTCbuoy(i,2)=max(meas);
        MTTbuoy(i,1)=min(simul);MTTbuoy(i,2)=min(meas);
        
        
        
        if GUIpp.PP3_timeshift==1
            if GUIpp.PP3_timeshiftBest==1
                [Xcorel,lags]=xcorr(meas,simul,'coeff');
               %[Xcorel,lags]=xcorr(simul,meas,'coeff');
                indMaxCorel=closest(Xcorel,max(Xcorel));
                OptCorr(i,2)=lags(indMaxCorel);
                %meas = circshift(meas,OptCorr(i,2));
                simul = circshift(simul,OptCorr(i,2));
             %   OptCorr(i,1)=(simul'*meas)/sqrt(varmeas*varsimul)/(length(meas));%correl(simul,meas);
                OptCorr(i,1)=(simul'*meas)/norm(simul)/norm(meas);
            else
                %meas = circshift(meas,GUIpp.PP3_timeshiftCustom);
                simul = circshift(simul,GUIpp.PP3_timeshiftCustom);
                %Corr(i)   = (simul'*meas)/sqrt(varmeas*varsimul)/(length(meas));%correl(simul,meas);
                Corr(i)   = (simul'*meas)/norm(simul)/norm(meas);
            end
          
        end
        
        simulSaved(:,i)=simul;
        measSaved(:,i)=meas;
        
        if any(bathy>0)
            simul=simul-mean(simul);
            meas=meas-mean(meas);
        end
        
        
        if GUIpp.PP3_Signals==1
            subaxis(Nfig,1,Nfig-i+1, 'Spacing', 0, 'Padding', 0, 'Margin', 0.13);
            
            plot(timeNow, meas,'b',timeNow, simul,'--r');
            
            
            if GUIpp.PPsetting_xlim==1
                xlim([GUIpp.PPsetting_xlimVal(1) GUIpp.PPsetting_xlimVal(2)])
            else
                xlim([time(indti) time(indtf)])
            end
            if GUIpp.PPsetting_ylim==1
                ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
            else
                ylim([-MM MM]);
            end
            
            
            %                         if mod(i,2)==0
            %                             set(gca,'ytick',[]);
            %                         else
            try
                if MM1*10>10
                    set(gca,'ytick',[-roundn(MM1,0) 0 roundn(MM1,0)]);
                elseif MM1*10>1
                    set(gca,'ytick',[-roundn(MM1,-1) 0 roundn(MM1,-1)]);
                else
                    set(gca,'ytick',[-roundn(MM1,-2) 0 roundn(MM1,-2)]);
                end
            catch
                
            end
            %                       end
            
            if i==1
                set(gca, 'box','off');
                xlabel('time [s]');
                %  legend('Meas','Simul')
            else
                set(gca,'xtick',[],'xcolor','w');
                set(gca,'box','off');
            end
            
            
            
            
            if i==round(Nfig/2) &&  mod(i,2)~=0
                if any(bathy>0)
                    ylabel('$\eta-\bar{\eta}$ [m]','Interpreter','latex');
                else
                    ylabel('\eta [m]');
                end
                FlagYlabel=1;
            end
            if FlagYlabel==0 && i==round(Nfig/2)+1
                if any(bathy>0)
                    ylabel('$\eta-\bar{\eta}$ [m]','Interpreter','latex');
                else
                    ylabel('\eta [m]');
                end
            end
            
            %th1=text(time(indti)+2, MM/2,['x=',num2str(roundn(MeasPos(i),-1))]);
            th1=text(time(indti)+1, MM/2,['W',num2str(i+1)]);
            plot_properties;
            
        end
        
    end

    
    if GUIpp.PP3_Signals==1
        set(hf3a, 'Position', [20 20 600 900])
        if GUIpp.save_figure==1
            if strcmp(GUIpp.save_figure_type,'.eps')
                saveas(hf3a,[sf_savename,'PP3_SimulMeasSignal',GUIpp.save_figure_type],'epsc')
            else
                saveas(hf3a,[sf_savename,'PP3_SimulMeasSignal',GUIpp.save_figure_type])
            end
        end
    end
    assignin('base','Corr',Corr)
    assignin('base','CorrN',CorrN)
    
    if  GUIpp.PP3_AsSk==1 ||...
            GUIpp.PP3_MTAMWL ==1
        if GUIpp.PPsetting_xlim==1
            indx1=closest(x,GUIpp.PPsetting_xlimVal(1));
            indx2=closest(x,GUIpp.PPsetting_xlimVal(2));
        else
            indx1=closest(x,xref(1)+Fbdy(1));
            indx2=closest(x,xref(2)-Fbdy(2));
        end
        indt1=closest(time,GUIpp.PP3_T(1));
        indt2=closest(time,GUIpp.PP3_T(2));
        
        if GUIpp.PP3_MTAMWL==1
            MWL_simul=mean(eta(indt1:dplot:indt2,indx1:indx2));
            MTC_simul=max(eta(indt1:dplot:indt2,indx1:indx2));
            MTT_simul=min(eta(indt1:dplot:indt2,indx1:indx2));
            
            hf=figure;
            set(hf, 'units','normalized','Position', [0.2, 0.2, 0.4,0.4]);
            plot(x(indx1:indx2),MTC_simul-MWL_simul,'-.k',x(indx1:indx2),MTT_simul-MWL_simul,':c',...
                x(indx1:indx2),MWL_simul,'--r',MeasPos,MTCbuoy(:,2)-MWL(:,2),'^b',...
                MeasPos,MTTbuoy(:,2)-MWL(:,2),'vb', MeasPos,MWL(:,2),'ob');
            legend('MTC-MWL (simul)','MTT-MWL ','MWL ','MTC-MWL (meas)','MTT-MWL ','MWL')
            xlim([x(indx1) x(indx2)])
            if GUIpp.PPsetting_ylim==1
                ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
            end
            xlabel('x[m]');ylabel('elevation [m]')
            plot_properties;
                    
            if GUIpp.save_figure==1
                if strcmp(GUIpp.save_figure_type,'.eps')
                    saveas(hf,[sf_savename,'PP3_MTA_MWL',GUIpp.save_figure_type],'epsc')
                else
                    saveas(hf,[sf_savename,'PP3_MTA_MWL',GUIpp.save_figure_type])
                end
            end
            
        end
        
        
   
        if GUIpp.PP3_AsSk==1
            Asym_simul=zeros(1,length(indx1:indx2));
            Skew_simul=Asym_simul;Hs_simul=Asym_simul;Kurt_simul=Hs_simul;
            
            [jbStop]=Java_stopbutton(statusbarObj);
            jProgressBar.setStringPainted( true );
            
            
            iter=1;Niter=indx2-indx1;tic;
            set(jProgressBar,'Maximum',Niter, 'Value',iter);
            statusbarObj.setText(['estimating time remaining..']);
            
            for i=indx1:indx2
                IDstop=eventLoopStop(jbStop);
                if IDstop==1, break;end;
                
                MeanSimul=mean(eta(indt1:dplot:indt2,i));
                Simul=eta(indt1:indt2,i);
                Hs_simul(i-indx1+1)= 4*sqrt(var(Simul,1));
                Hilb_sim=imag(hilbert(Simul-MeanSimul));
                Asym_simul(i-indx1+1)=mean(Hilb_sim.^3)/((mean((Simul-MeanSimul).^2))^(3/2));
                Skew_simul(i-indx1+1)=mean((Simul-MeanSimul).^3)/((mean((Simul-MeanSimul).^2))^(3/2));
                Kurt_simul(i-indx1+1)=kurtosis(Simul);
                
                if mod(iter,floor(0.1*Niter))==0 || iter==floor(0.01*Niter)
                    set(jProgressBar,'Maximum',Niter, 'Value',iter);
                    ETA=remain_time(iter,Niter);
                    statusbarObj.setText(['time remaining=', num2str(ETA)]);
                end
                iter=iter+1;
            end
            statusbarObj.setText(['']);
            
            
            
            hf=figure;
            set(hf, 'units','normalized','Position', [0.2, 0.2, 0.4,0.4]);
            subplot(4,1,1)
            plot(x(indx1:stepdx:indx2),Hs_simul(1:stepdx:end),'r',MeasPos,Hs(:,2),'ob')
            xlabel('x[m]');
            ylabel('Hs [m]')
            plot_properties;
            xlim([x(indx1) x(indx2)])
            if GUIpp.PPsetting_ylim==1
                ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
            end
            subplot(4,1,2)
            plot(x(indx1:stepdx:indx2),Skew_simul(1:stepdx:end),'r',MeasPos,Skew(:,2),'ob')
            xlabel('x[m]');
            ylabel('Sk ')
            plot_properties;
            xlim([x(indx1) x(indx2)])
            if GUIpp.PPsetting_ylim==1
                ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
            end
            subplot(4,1,3)
            plot(x(indx1:stepdx:indx2),-Asym_simul(1:stepdx:end),'r',MeasPos,-Asym(:,2),'ob')
            xlabel('x[m]');
            ylabel('-Asym')
            plot_properties;
            xlim([x(indx1) x(indx2)])
            if GUIpp.PPsetting_ylim==1
                ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
            end
            subplot(4,1,4)
            plot(x(indx1:stepdx:indx2),Kurt_simul(1:stepdx:end),'r',MeasPos,kurt(:,2),'ob')
            xlabel('x[m]');
            ylabel('Kurt')
            plot_properties;
            xlim([x(indx1) x(indx2)])
            if GUIpp.PPsetting_ylim==1
                ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
            end
            
            
            if GUIpp.save_figure==1
                if strcmp(GUIpp.save_figure_type,'.eps')
                    saveas(hf,[sf_savename,'PP3_Hs_Sk_Asym',GUIpp.save_figure_type],'epsc')
                else
                    saveas(hf,[sf_savename,'PP3_Hs_Sk_Asym',GUIpp.save_figure_type])
                end
            end
            
            
            set(jProgressBar,'Maximum',Niter, 'Value',Niter);
            ETA=remain_time(Niter,Niter);
            statusbarObj.setText(['time remaining=', num2str(ETA)]);
            jProgressBar.setVisible(0);
            jbStop.setVisible(0);
            statusbarObj.setText('done.');
            
        end
    end
    
    
    quantData.vari=vari; quantData.Corr=Corr; quantData.kurt=kurt;
    quantData.Hs=Hs;   quantData.OptCorr=OptCorr;
    quantData.Asym=Asym; quantData.Skew=Skew;
    
    if GUIpp.PP3_AsSk==1
        quantData_xinterv.Asym=Asym_simul;
        quantData_xinterv.Skew=Skew_simul;
        quantData_xinterv.Hs=Hs_simul;
        quantData_xinterv.Kurt=Kurt_simul;
    else
        quantData_xinterv=[];
    end
    
    PPdata.PP3.quantDataSimul_xinterv=quantData_xinterv;
    PPdata.PP3.quantDataAtbuoy=quantData;
    PPdata.PP3.simul=simulSaved;
    PPdata.PP3.meas=measSaved;
    PPdata.PP3.time=timeNow;
    PPdata.PP3.MeasPosition=MeasPos;
    PPdata.PP3.spSIMUL=spSIMUL;
    PPdata.PP3.spMEAS=spMEAS;
    PPdata.PP3.frequency_rad=ww;
    
    
    assignin('base','quantData',quantData);
    if GUIpp.PP3_timeshift==1
        if GUIpp.PP3_timeshiftBest==1
            XCORR=OptCorr(:,1);
        else
            XCORR=Corr;
        end
    else
        XCORR=Corr;
    end
    
    
    
    if GUIpp.PP3_Quant==1
        %f = figure('Position',[440 500 461 146]);
        f = figure('units','normalized','Position',[0.65 0.5 0.25 0.13],'name','Quantitative information.');
        
        var3=vari(:,3);kurt3=kurt(:,3);Hs3=Hs(:,3);Sk3=Skew(:,3);As3=Asym(:,3);
        d = [XCORR';var3';Hs3';kurt3';Sk3';As3'];
        
        %        Create the column and row names in cell arrays
        for i=1:Nfig
            cnames{i} = ['W_',num2str(i+1)];
        end
        
        rnames = {'Corr','Var','Hs','Kurt','Sk','As'};
        %      Create the uitable
        t = uitable(f,'units','normalized','Position',[0 0 1 1],'Data',d,...
            'ColumnName',cnames,...
            'RowName',rnames,'ColumnWidth','auto');
        title('Table. Quantitative information [Simul/Meas]')
        if GUIpp.save_figure==1
            saveas(f,[sf_savename,'PP3_quantitative','.fig'])
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%Plot Spectra%%%%%%%%%%%%%%%
    Nw=length(ww);indEn=floor(Nw/2); %6 minimal cutfracwn
    if ESpId~=2
        [~,sppeakindexPP]=max(abs(spSIMUL(:,1)));
        [MMsp,~]=max(max(spSIMUL));
    else
        [MMsp,sppeakindexPP]=max(abs(spSIMUL(:,1)));
        [MMspEn,sppeakindexPPEn]=max(abs(spSIMULEn(:,1)));
        if MMspEn>MMsp, MMsp=MMspEn; sppeakindexPP=sppeakindexPPEn;end;
    end
    
    if GUIpp.PP3_Spectra==1
        hf3b=figure;
        FlagYlabel=0;
        for i=1:Nfig
            subaxis(Nfig,1,Nfig-i+1, 'Spacing', 0, 'Padding', 0, 'Margin', 0.13);
            
            if ESpId~=2
                plot(ww(1:dplot:indEn), spMEAS(1:dplot:indEn,i),'b',ww(1:dplot:indEn), spSIMUL(1:dplot:indEn,i),'--r');
                
                if i==1
                    set(gca, 'box','off');
                    xlabel('\omega [rad/s]');
                else
                    set(gca,'xtick',[],'xcolor','w');
                    set(gca,'box','off');
                end
            else
                plot(ww(1:dplot:indEn), spMEAS(1:dplot:indEn,i),'b',ww(1:dplot:indEn), spMEASEn(1:dplot:indEn,i),'--b',...
                    ww(1:dplot:indEn), spSIMUL(1:dplot:indEn,i),'r',ww(1:dplot:indEn), spSIMULEn(1:dplot:indEn,i),'--r');
                
                if i==1
                    set(gca, 'box','off');
                    xlabel('\omega [rad/s]');
                    legend('Meas','Simul')
                else
                    set(gca,'xtick',[],'xcolor','w');
                    set(gca,'box','off');
                end
            end
            
            nupeakPP=max(influx.nu_p,ww(sppeakindexPP));
            
            
            
            if GUIpp.PPsetting_xlim==1
                xlim([GUIpp.PPsetting_xlimVal(1) GUIpp.PPsetting_xlimVal(2)])
            else
                xlim([0 min(nupeakPP*7,ww(indEn))]);
            end
            if GUIpp.PPsetting_ylim==1
                ylim([GUIpp.PPsetting_ylimVal(1) GUIpp.PPsetting_ylimVal(2)])
            else
                ylim([0 MMsp+MMsp/4]);
            end
            
            
            try
                %                 if mod(i,2)==0
                %                     set(gca,'ytick',[]);
                %                 else
                if MMsp*10>10
                    Y1=roundn(MMsp/2,0);Y2=roundn(MMsp,0);
                    if Y1~=Y2
                        set(gca,'YTick',[0 roundn(MMsp/2,0) roundn(MMsp,0) ])
                    else
                        set(gca,'YTick',[0  roundn(MMsp,0) ])
                    end
                elseif MMsp*10>1
                    Y1=roundn(MMsp/2,-1);Y2=roundn(MMsp,-1);
                    if Y1~=Y2
                        set(gca,'YTick',[0 roundn(MMsp/2,-1) roundn(MMsp,-1) ])
                    else
                        set(gca,'YTick',[0  roundn(MMsp,-1) ])
                    end
                else
                    Y1=roundn(MMsp/2,-2);Y2=roundn(MMsp,-2);
                    if Y1~=Y2
                        set(gca,'YTick',[0 roundn(MMsp/2,-2) roundn(MMsp,-2) ])
                    else
                        set(gca,'YTick',[0  roundn(MMsp,-2) ])
                    end
                end
                %end
            catch
                
            end
            
            
            if i==round(Nfig/2) &&  mod(i,2)~=0
                if ESpId==0,ylabel('Amplitude Sp.');
                else ylabel('Variance dens. Sp.');end;
                FlagYlabel=1;
            end
            if FlagYlabel==0 && i==round(Nfig/2)+1
                if ESpId==0,ylabel('Amplitude Sp.');
                else ylabel('Variance dens. Sp.');end;
            end
            
            % th1=text(ww(5), roundn(spSIMUL(sppeakindexPP,1)/2,-2),['x=',num2str(roundn(MeasPos(i),-1))]);
            th1=text(ww(5), roundn(spSIMUL(sppeakindexPP,1)/2,-2),['W',num2str(i+1)]);
            
            plot_properties;
            %         set(p0, 'linewidth', 0.5);
        end
        set(hf3b, 'Position', [600 20 600 900])
        if GUIpp.save_figure==1
            if strcmp(GUIpp.save_figure_type,'.eps')
                saveas(hf3b,[sf_savename,'PP3_SimulMeasSpectra',GUIpp.save_figure_type],'epsc');
            else
                saveas(hf3b,[sf_savename,'PP3_SimulMeasSpectra',GUIpp.save_figure_type]);
            end
        end
    end
end

if GUIpp.PPsetting_save_data==1
    save([sf_savename,'PPdata'], 'PPdata');
end

% temp=zeros(length(simulSaved(:,1))+1,length(simulSaved(1,:))+1);
% temp(1,2:end)=MeasPos;
% temp(2:end,1)=timeNow;
% temp(2:end,2:end)=simulSaved;
% save([sf_savename,GUIpp.savename,'_simul_compare.txt'], 'temp','-ascii');
