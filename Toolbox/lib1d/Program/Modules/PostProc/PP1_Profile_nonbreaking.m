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

if shippar.check==1
    IndSL=closest(x,shippar.form.xShip(1));
    IndSR=closest(x,shippar.form.xShip(3));
    
    MTT(IndSL+1:IndSR-1)=NaN;
    MTC(IndSL+1:IndSR-1)=NaN;
    %             figure;
    %             plot(x,MTT,'r',x,MTC,'b')
    %
    draft=max(shippar.form.Sdraft);
    slength=max(shippar.form.Slength);
   
    z0SS=linspace(-draft,1.2*max(MTC),20);
    x0SS=linspace(-slength/2,slength/2,20);
    onez=ones(size(z0SS));
    onex=ones(size(x0SS));
    
    zSS=[z0SS       -draft.*onex    z0SS            1.2*max(MTC).*onex];
    xSS=[x0SS(1)*onez x0SS          x0SS(end).*onez fliplr(x0SS)];
    xS0=shippar.form.xShip(2);
    
    
    
    %             figure;
    %             fill(xSSn0,zSSn0,'r')
end



for I=1:length(tsnap)
         if idCol==9, idCol=1; end;
         indt=(closest(time,tsnap(I)));
         plot(x(1:dplot:end), eta(indt,1:dplot:end),CR{idCol});
         idCol=idCol+1;
         leg{I}=['\eta t=',num2str(tsnap(I)),'[s]'];
         hold on;
         
         if shippar.check==1
             sZ=shipRB(indt,3);
             sX=shipRB(indt,1);
             sTheta=shipRB(indt,5);
             xSSn=xS0+sX+(xSS).*cos(sTheta)-(zSS).*sin(sTheta);
             zSSn=sZ+(xSS).*sin(sTheta)+(zSS).*cos(sTheta);
             
             fill(xSSn,zSSn,'b')
           
         end
         
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

if IdOnlyEta==1&&length(tsnap)==1
else
%legend(leg{1:end},'Location','northwest')
end
%title(title_text);

PPdata.PP1.Profile_x=x(1:dplot:end);
PPdata.PP1.Profile_eta=eta(closest(time,tsnap),1:dplot:end);
PPdata.PP1.Profile_tsnap=tsnap;
PPdata.PP1.Profile_MTC=MTC;
PPdata.PP1.Profile_MTT=MTT;

