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
%%%%%%%%%%%%%%%%%  Pre-Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
          
    SpatialSetup;          % Make spatial interval, incl bathy, bdy, wall
    GenerationSetup;       % Prepare influx signals for embedded influx source
    if shippar.check==1
    ShipSetup;
    end
    OperatorSetup;
    if shippar.check==1
    ShipPotentialSetup;
    end
    
    parambreak;
   
   
   
    
    PrepView;              % Preview
    ODE_partition_setup;
     
    if strcmp(options.log,'Yes')
    logfile_display;
    end  

          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%Save to handle of preproc pushbutton%%%%%
preproc.influx=influx;   preproc.par     =par;
preproc.model =model;    preproc.Proj    =Proj;
preproc.bath  =bath;     preproc.input   =input;
preproc.options=options; preproc.IVP     =IVP;
preproc.Oprt   =Oprt;    preproc.shippar =shippar;
preproc.bdyassim=bdyassim;
if options.IdFlagGui==1
set(handles.pushbutton1,'userdata',preproc);    
set(jProgressBar,'Maximum',100, 'Value',100);
jProgressBar.setStringPainted( true );    
end

if ~isdeployed
assignin('base','Proj',Proj)
assignin('base','IVP',IVP)
assignin('base','par',par)
assignin('base','influx',influx)
assignin('base','model',model)
assignin('base','options',options)
assignin('base','input',input)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
