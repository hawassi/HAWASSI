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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%    HAWASSI 1D                                          %%%%%%%%%%
%%%%%%%%%    Hamiltonian Wave-Ship-Structure Interaction         %%%%%%%%%%
%%%%%%%%%    copyright (c): LabMath-Indonesia & Unversity Twente %%%%%%%%%%
%%%%%%%%%    version 23 June 2015                                %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  ODE-Solve %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tintval      =influx.gen.timesig';

 
%%%%%%%%%%%%%%%%%%%%%%%%%%Evolution%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
evol_partition;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
par.x      =x;
par.CPUtime=Timeodesolver;
par.Crel=CompRel;
log_txt.par=par;

if ~isdeployed
assignin('base','Proj',Proj)
assignin('base','bath',bath)
assignin('base','par',par)   %copy to workspace
assignin('base','options',options)
assignin('base','influx',influx)
assignin('base','CompRel',CompRel);  %copy to workspace
assignin('base','Timeodesolver',Timeodesolver);
end

if strcmp(options.log,'Yes')
    diary on;
    disp('After simulation:');
    disp(['CPUtime for ODEs: ',num2str(Timeodesolver,'%10.2e\n'),' sec']);
    disp(['Relative computation time: ',num2str(CompRel,2),' percent']);
    disp(['Data saved as : ',Proj.savename,'_simul.mat'])
end
comptime.cpu=Timeodesolver;
comptime.crel=CompRel;
if options.IdFlagGui==1
set(handles.togglebutton_logfile,'userdata',comptime);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%Calling Post-Proc GUI and/or Interior Flow GUI%%%%%%%%%%%%%%%%%
if options.IdFlagGui==1
    statusbarObj.setText('Post-processing');
    PostProc_GUI(inputPProc)
    
    
    if options.interior.check==1
        if checkInterior==1
            InteriorProp_GUI(InputInteriorCalc);
        end
    end
end
global Idstop 
if time(end)<tintval(end)
 if Idstop==0
 statusbarObj.setText('time-stepping failure!');
 disp('---------------------------------------------------------------------');
 disp('Warning: time-stepping failure')
 disp('---------------------------------------------------------------------');
 statusbarTxt = statusbarObj.getComponent(0);
 statusbarTxt.setForeground(java.awt.Color.red);    
 else
 statusbarObj.setText('terminated by user.');
 statusbarTxt = statusbarObj.getComponent(0);
 statusbarTxt.setForeground(java.awt.Color.blue);    
 end
else
    statusbarObj.setText('done.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.blue);
end

if input.assim.check==0
    if options.partition.combine==1 && Npartition>1
        if Iddelete==1
            statusbarObj.setText('done.');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.blue);
        else
            statusbarObj.setText('Combining file failed. Data is too large');
            disp('---------------------------------------------------------------------');
            disp('Combining file failed. Data is too large')
            disp('---------------------------------------------------------------------');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
        end
    end
else
     if options.partition.combine==1 && ITERAssimSave-1>1
        if Iddelete==1
            statusbarObj.setText('done.');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.blue);
        else
            statusbarObj.setText('Combining file failed. Data is too large');
            disp('---------------------------------------------------------------------');
            disp('Combining file failed. Data is too large')
            disp('---------------------------------------------------------------------');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
        end
    end
end
if Flag_timestep_failure==1
    statusbarObj.setText('time-stepping failure!');
    disp('---------------------------------------------------------------------');
    disp('Warning: time-stepping failure')
    disp('---------------------------------------------------------------------');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
end
if strcmp(options.log,'Yes')
diary off;
end
if options.IdFlagGui==1
    if shippar.check==1
        set(handles.pushbutton_ship_setupGui,'Enable','on')
    end
    clearvars -except handles;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Done%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%