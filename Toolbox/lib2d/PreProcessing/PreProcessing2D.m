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
%%%%%%%%%% HAWASSI-AB 2D                                         %%%%%%%%%%
%%%%%%%%%% Hamiltonian Wave-Ship-Structure Interaction           %%%%%%%%%%
%%%%%%%%%% Copyright (c): LabMath-Indonesia                      %%%%%%%%%%
%%%%%%%%%% version: 5 July 2016                                  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RK%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Pre-processing                                        %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.mc.check==0
    set(jProgressBar,'Maximum',100, 'Value',10);
    SpatialSetup2D;          % Make spatial interval, incl bathy, bdy, wall
    set(jProgressBar,'Maximum',100, 'Value',30);
    GenerationSetup2D;       % Prepare influx signals for embedded influx source
    set(jProgressBar,'Maximum',100, 'Value',60);
    OperatorSetup2D;         % Preparing Operators that will be applied in RHS
    ParamBreakSetup2D;       % Preparing parameters for breaking model
    %PrepView;              % Preview
    ODE_partition_setup2D;
    set(jProgressBar,'Maximum',100, 'Value',80);
    if strcmp(options.log,'Yes')
        logfile_display2D;
    end
    
    %%%%%%%%%%%%%%%%%Passing variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Preproc.log=log;
    Preproc.par=par;Preproc.model=model;
    Preproc.bath=bath;Preproc.input=input;
    Preproc.Proj=Proj;Preproc.timeSimul=timeSimul;
    Preproc.dom=dom;Preproc.options=options;
    Preproc.spatial=spatial;
    Preproc.ivp=ivp;Preproc.Oprt=Oprt;
    Preproc.parBreak=parBreak;Preproc.influx=influx;
    Preproc.body=input.body;
    Preproc.bdyassim=bdyassim;
else
    Preproc(1:options.mc.N+1)=struct;
    propdata=input.wave.propertiesdata;
    mcwaveinput=options.mc.waveinput_data;
    Ninf=length(propdata(:,1));
    for JJ=1:options.mc.N+1
        if JJ>1
            if options.mc.numset_check==1
                Npx=length(options.mc.numset_px(:));
                if JJ>Npx
                    spatial.dom.dx=options.mc.numset_px(Npx);
                else
                    spatial.dom.dx=options.mc.numset_px(JJ-1);
                end
                Npy=length(options.mc.numset_py(:));
                if JJ>Npy
                    spatial.dom.dy=options.mc.numset_py(Npy);
                else
                    spatial.dom.dy=options.mc.numset_p(JJ-1);
                end
            end
        end    
            if options.mc.waveinput_check==1
                if JJ>1
                    for ii=1:Ninf
                        if strcmpi(propdata(ii,1),'Harmonic')
                            if isfield(mcwaveinput(ii),'A_param')
                                NN=length(mcwaveinput(ii).A_param(:));
                                if JJ>NN
                                propdata(ii,2)={mcwaveinput(ii).A_param(end)};
                                else    
                                propdata(ii,2)={mcwaveinput(ii).A_param(JJ-1)};
                                end
                            end
                            if isfield(mcwaveinput(ii),'Tp_param')
                                NN=length(mcwaveinput(ii).Tp_param(:));
                                if JJ>NN
                                propdata(ii,4)={mcwaveinput(ii).Tp_param(end)};
                                else
                                propdata(ii,4)={mcwaveinput(ii).Tp_param(JJ-1)};
                                end
                            end
                        elseif strcmpi(propdata(ii,1),'Jonswap')
                            if isfield(mcwaveinput(ii),'Hs_param')
                                NN=length(mcwaveinput(ii).Hs_param(:));
                                if JJ>NN
                                 propdata(ii,3)={mcwaveinput(ii).Hs_param(end)};
                                else
                                propdata(ii,3)={mcwaveinput(ii).Hs_param(JJ-1)};
                                end
                            end
                            if isfield(mcwaveinput(ii),'Tp_param')
                                NN=length(mcwaveinput(ii).Tp_param(:));
                                if JJ>NN
                                propdata(ii,4)={mcwaveinput(ii).Tp_param(end)};
                                else
                                propdata(ii,4)={mcwaveinput(ii).Tp_param(JJ-1)};
                                end
                            end
                            if isfield(mcwaveinput(ii),'gamma_param')
                                NN=length(mcwaveinput(ii).gamma_param(:));
                                if JJ>NN
                                propdata(ii,5)={mcwaveinput(ii).gamma_param(end)};
                                else
                                propdata(ii,5)={mcwaveinput(ii).gamma_param(JJ-1)};
                                end
                            end
                            if isfield(mcwaveinput(ii),'s_param')
                                NN=length(mcwaveinput(ii).s_param(:));
                                if JJ>NN
                                propdata(ii,6)={mcwaveinput(ii).s_param(end)};
                                else
                                propdata(ii,6)={mcwaveinput(ii).s_param(JJ-1)};
                                end
                            end
                        elseif strcmpi(propdata(ii,1),'User defined (variance density spectrum)')
                            if isfield(mcwaveinput(ii),'s_param')
                                NN=length(mcwaveinput(ii).s_param(:));
                                if JJ>NN
                                propdata(ii,6)={mcwaveinput(ii).s_param(end)};
                                else
                                propdata(ii,6)={mcwaveinput(ii).s_param(JJ-1)};
                                end
                            end
                        end
                    end
                end
                input.wave.propertiesdata=propdata;
            end
            
            set(jProgressBar,'Maximum',100, 'Value',10);
            try
                SpatialSetup2D;          % Make spatial interval, incl bathy, bdy, wall
            catch
                 statusbarObj.setText('Error in the SpatialSetup2D.'); 
            end
            set(jProgressBar,'Maximum',100, 'Value',30);
            try
            GenerationSetup2D;       % Prepare influx signals for embedded influx source
            catch
            statusbarObj.setText('Error in the GenerationSetup2D.');    
            end
            set(jProgressBar,'Maximum',100, 'Value',60);
            try
            OperatorSetup2D;         % Preparing Operators that will be applied in RHS
            catch
            statusbarObj.setText('Error in the OperatorSetup2D.');       
            end
            try
            ParamBreakSetup2D;       % Preparing parameters for breaking model
            catch
            statusbarObj.setText('Error in the ParamBreakSetup2D.');    
            end
           
            ODE_partition_setup2D;
            set(jProgressBar,'Maximum',100, 'Value',80);
            if strcmp(options.log,'Yes')
                logfile_display2D;
            end
        
            %%%%%%%%%%%%%%%%%Passing variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Preproc(JJ).log=log;
            Preproc(JJ).par=par;
            Preproc(JJ).model=model;
            Preproc(JJ).bath=bath;
            Preproc(JJ).input=input;
            Preproc(JJ).Proj=Proj;
            Preproc(JJ).timeSimul=timeSimul;
            Preproc(JJ).dom=dom;
            Preproc(JJ).options=options;
            Preproc(JJ).spatial=spatial;
            Preproc(JJ).ivp=ivp;
            Preproc(JJ).bdyassim=bdyassim;
            Preproc(JJ).Oprt=Oprt;
            Preproc(JJ).parBreak=parBreak;
            Preproc(JJ).influx=influx;
            Preproc(JJ).body=input.body;
        
    end      
end