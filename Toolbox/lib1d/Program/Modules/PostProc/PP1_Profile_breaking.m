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
IdOnlyEta=0;
for I=1:length(tsnap)
    if idCol==9, idCol=1; end;   
    plot(x(1:dplot:end), eta(closest(time,tsnap(I)),1:dplot:end),CR{idCol});
    idCol=idCol+1;
    leg{I}=['\eta t=',num2str(tsnap(I)),'[s]'];
    hold on;
end
if I==1
leg{I}=['\eta'];    
end
IdLeg=I;       
if GUIpp.PP1_Prof_MTA==1 && GUIpp.PP1_Prof_Bathy==0
        plot(x(1:dplot:end), MTC(1:dplot:end), '-.k',...
                x(1:dplot:end),MTT(1:dplot:end), ':c');
        leg{IdLeg+1}='MTC';leg{IdLeg+2}='MTT';
        Nt=length(tsnap);
        if Nt~=0
            if Nt~=1
                for i=1:Nt-1
                    tsnap_str{i}=[num2str(tsnap(i)),';',];
                end
            end
            tsnap_str{Nt}=num2str(tsnap(Nt));
        end

        title_text=[dynmodel,', MTA & profile @ time: ',[tsnap_str{1:Nt}], ' [s]'];
        IdLeg=IdLeg+2;
elseif GUIpp.PP1_Prof_MTA==0 && GUIpp.PP1_Prof_Bathy==1
    scale=GUIpp.PP1_Prof_Bathy_scale;
    plot(x(1:dplot:end), bathy(1:dplot:end)*scale, '--k');
    leg{IdLeg+1}=['Bathy, scale: ',num2str(scale)];
    IdLeg=IdLeg+1;
    Nt=length(tsnap);
    if Nt~=0
        if Nt~=1
            for i=1:Nt-1
                tsnap_str{i}=[num2str(tsnap(i)),';',];
            end
        end
        tsnap_str{Nt}=num2str(tsnap(Nt));
    end

    title_text=[dynmodel,', Bathymetry & profile @ time: ',[tsnap_str{1:Nt}], ' [s]'];
elseif GUIpp.PP1_Prof_MTA==1  && GUIpp.PP1_Prof_Bathy==1
    plot(x(1:dplot:end), MTC(1:dplot:end), '-.k',...
                x(1:dplot:end),MTT(1:dplot:end), ':c');
        leg{IdLeg+1}='MTC';leg{IdLeg+2}='MTT';
    hold on;
    scale=GUIpp.PP1_Prof_Bathy_scale;
    plot(x(1:dplot:end), bathy(1:dplot:end)*scale, '--k');
    leg{IdLeg+3}=['Bathy, scale: ',num2str(scale)];
    IdLeg=IdLeg+3;
     Nt=length(tsnap);
    if Nt~=0
        if Nt~=1
            for i=1:Nt-1
                tsnap_str{i}=[num2str(tsnap(i)),';',];
            end
        end
        tsnap_str{Nt}=num2str(tsnap(Nt));
    end
    title_text=[dynmodel,', Bathymetry, MTA & profile @ time: ',[tsnap_str{1:Nt}], ' [s]'];
else
     Nt=length(tsnap);
    if Nt~=0
        if Nt~=1
            for i=1:Nt-1
                tsnap_str{i}=[num2str(tsnap(i)),';',];
            end
        end
        tsnap_str{Nt}=num2str(tsnap(Nt));
    end
    title_text=[dynmodel,', profile @ time ',[tsnap_str{1:Nt}],' [s]'];
    IdOnlyEta=1;
end

if strcmp(wall.presence,'Yes')
    hold on;
    plot(Xwall,(-2*ampli:ampli/5:2*ampli),'oy');
    leg{IdLeg+1}='Wall';IdLeg=IdLeg+1;
end

dt=time(2)-time(1);
    
for I=1:length(tsnap)    
if time(closest(time,tsnap(I)))>=dataBreak_nodes(1,1)&&...
        time(closest(time,tsnap(I)))<=dataBreak_nodes(end,1)
    indtBreak=closest(dataBreak_nodes(:,1),time(closest(time,tsnap(I))));
    indendb=find(dataBreak_nodes(indtBreak,:)==0,1,'first')-1;
    indxBreak=dataBreak_nodes(indtBreak,2:indendb);
    hold on;
    if all(abs(tsnap(I)-dataBreak_nodes(:,1))>dt)   
    else
    plot(x([indxBreak]),eta(closest(time,tsnap(I)),[indxBreak]),'ob');
    end
end 
end

if IdOnlyEta==1&&length(tsnap)==1
else
%legend(leg{1:end},'Location','northwest')
end
%title(title_text);


PPdata.PP1.Profile_x=x(1:dplot:end);
PPdata.PP1.Profile_eta=eta(closest(time,tsnap),1:dplot:end);
% PPdata.PP1.Profile_xbreak=x([indxBreak]);
% PPdata.PP1.Profile_etabreak=eta(closest(time,tsnap),[indxBreak]);
PPdata.PP1.Profile_tsnap=tsnap;
PPdata.PP1.Profile_MTC=MTC;
PPdata.PP1.Profile_MTT=MTT;
