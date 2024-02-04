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
%%%%%%%%%%%%%%%%%%%%%%%%%%%  MAIN-USER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%Project Preparation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Proj.path                = GUImain.proj.path;
Proj.projdir             = GUImain.proj.projdir;
Proj.Dir                 = GUImain.proj.workdir;  % working directory
Proj.savename            = GUImain.proj.name;    % save simulated data as 'savename'.mat
Proj.UseNote             = GUImain.proj.note;        % This UserNote is inserted in outputlog

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model.evol             = GUImain.evol;        % 'HS' for dispersive Ham. System
% 'SWE' for shallow water equation
model.nonlinear        =  GUImain.nonlin;      % for nonlin order: #=1(linear),2,3,4;
model.breaking.check   =  GUImain.Break;       % 'Yes' or 'No'
if strcmp(model.breaking.check,'Yes')
    model.breaking.KBC         = GUImain.Break_param(1);    % Kin.Br.Cond.: U/C in [0.7:1]
    if  GUImain.BreakDef==0
        model.breaking.TC      = GUImain.Break_paramDef(1); % Terminal Condition uF=TC*uI; TC in [0.2:0.3]
        model.breaking.Tchar   = GUImain.Break_paramDef(2); % Characteristic Time T*=Tchar*Tp; Tchar in [0.1:0.5]; Tp->peak period
    else  %default
        model.breaking.TC      =0.2;
        model.breaking.Tchar   =0.5;
    end
    model.breaking.delb    = 1.2;              % Mixing length Coef
end
model.dispersion       =  GUImain.disp;        % OmExact,OmSWE,OmKdV,OmBBM,  available; specify OmUser
model.groupvel         =  GUImain.groupvel;    % specify UgUser related to OmUser
model.detaC2            = GUImain.detaC2;       % dC2/deta for run-up case
model.Cder             = GUImain.Cder;        % dC/deta for run-up case
model.OmFun            =  GUImain.dispEqInput;
model.UgFun            =  GUImain.groupVelEqInput;

model.influx.type      =  GUImain.inftype ;    % 'Area', 'Area-short', 'Point'
model.influx.direction =  GUImain.infdir;      % = 'Uni+' for unidrectional to right
% = 'Uni-' for unidrectional to left
% = 'Bi' for symmetric bidirectional both sides
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%Wave-Input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%Initial Conditions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IVP.type               =GUImain.IVP.type;
IVP.typename           =GUImain.IVP.typename;
IVP.A                  =GUImain.IVP.A;
IVP.lambda             =GUImain.IVP.lambda;
IVP.x0                 =GUImain.IVP.x0;
if IVP.type==4
    IVP.data               =GUImain.IVP.file.data;
    IVP.filename           =GUImain.IVP.file.name;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%Influxing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input.type              =  GUImain.WT;
if input.type==2
    input.Tp            =  GUImain.Tp;         % [s] Peak period
    input.Hs            =  GUImain.Hs;         % [m] significant wave height
    input.name          =  'Harmonic';
elseif input.type==3
    input.Tp            =  GUImain.Tp;
    input.Hs            =  GUImain.Hs;
    input.JS_gamma      =  GUImain.Jsg;        %  parameter gamma in JS-spectrum
    input.name          =  'Jonswap';
elseif input.type==4
    input.usersignal    = GUImain.WT_O_data(1:end,2); %user data
    input.timesig       = GUImain.WT_O_data(1:end,1);
    input.name          = 'User-defined';
else
    input.name          = 'None';
end
input.filter.check=GUImain.Filter;
if input.filter.check==1
    input.filter.LFreq=GUImain.Filter_LFHF(1); %Low frequency filter
    input.filter.HFreq=GUImain.Filter_LFHF(2); %High frequency filter
end
input.ramp.check      =  GUImain.ramp.check;   %
input.ramp.length     =  GUImain.ramp.val;     % the length of the ramp function #Tp

%%%%%%%%%%%%%%%%%%%%%%%%Assimilation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input.assim.check     =GUImain.Assim.check;
input.assim.data      =GUImain.Assim.data;
%%%%%%%%%%%%%%%%%%%%%%%%simulation time%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

input.dt              =  GUImain.tstep;        % time step dt [s] to assemble ODE output
input.t_init          =  GUImain.tinterv(1);   % simulation start time [s]
input.t_end           =  GUImain.tinterv(2);   % simulation end time [s]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%Wind-input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input.wind.check      = GUImain.wind.check;
if input.wind.check==1
    input.wind.coef      = GUImain.wind.coef;
    input.wind.tinterv   = GUImain.wind.tinterv;
    input.wind.xinterv   = GUImain.wind.xinterv;
else
    input.wind.coef      = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%Spatial-Space%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spatial.Xleft           = GUImain.xinterv(1);  % Start position
spatial.Xright          = GUImain.xinterv(2);  % End   position
spatial.Xinflux         = GUImain.Xinflux;     % influx position
spatial.dx              = GUImain.dx;           % # spatial discritization 2^p;
spatial.cutfrac         = GUImain.cutfracK;    % anti-aliassing fraction factor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%Bathymetry%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bath.type=GUImain.bathy;
if strcmp(bath.type,'F')==1
    bath.depth          = GUImain.depth;       %depth in [m]
    bath.name           = 'Flat';
elseif strcmp(bath.type,'B')==1
    bath.par            = [GUImain.depth(1);GUImain.depth(2);GUImain.slope(1);GUImain.slope(2)];%[30,10,1/10,200]; % for slope: par.bath = [Depth_left,Depth_right,slope,X_startslope];
    bath.name           = 'Slope';
    bath.interpRef      = GUImain.bathyInterp+1;
    bath.DmidRef        = GUImain.DmidRef;
elseif strcmp(bath.type,'BR')==1          %bathy runup
    bath.par            = [GUImain.depth(1);GUImain.slope(1);GUImain.slope(2)];
    bath.hmin           = GUImain.depth(2); % for runup
    bath.name           = 'Slope (Shore)';
    bath.interpRef      = GUImain.bathyInterp+1;
    bath.DmidRef        = GUImain.DmidRef;
else
    bath.name           = 'User-defined';
    bath.data           = GUImain.bathyO.data;% 'Bathymetry data'
    bath.filename       = GUImain.bathyO.filename;% 'Bathymetry data'
    bath.interpRef      = GUImain.bathyInterp+1;
    bath.DmidRef        = GUImain.DmidRef;
    if any(bath.data>=0)
        bath.hmin           = GUImain.depth(1); % for runup
    end
    
end
bath.friction.check     = GUImain.bottomfriction_check;
bath.friction.data      = GUImain.frictiondata;


%%%WALL
bath.wall.check         = GUImain.Wall;         % No or Yes
if strcmp(bath.wall.check,'Yes')
    bath.wall.length        = 20;%
    bath.wall.position      = GUImain.Wall_position; % [m]
    bath.wall.method        = GUImain.Wall_method;
    bath.wall.type          = GUImain.Wall_type;
    % = 1: reflPercent to be specified
    % = 2: user-defined reflection per frequency
    bath.wall.Coef           = GUImain.Wall_Coef;
    bath.wall.file_R         = GUImain.Wall_reflEq;
    bath.wall.Evmodes        = 2;
end

bath.FBL                 = GUImain.FBL;          % [m] length interval at both ends to prevent Fourier looping/ discontinuities
bath.influx_AdjZone      = GUImain.NonlinAdj;   % length (in # of  peak-wavelength) for adjustment-zone
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% Ship Setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shippar                  = GUImain.shippar;      %%shippar.check shippar.data shippar.user_shapedata
if shippar.check==1
    shippar.Evmode           = shippar.Evmode;
    shippar.moored.check     = cell2mat(shippar.data(1,11));
    shippar.calclin          = shippar.linearCal;   
    try
        shippar.moored.Tn1        = str2num(cell2mat(shippar.data(1,12)));%
        shippar.moored.Tn3        = str2num(cell2mat(shippar.data(1,13)));%
        shippar.moored.Tn5        = str2num(cell2mat(shippar.data(1,14)));%
    catch
        shippar.moored.Tn1        = cell2mat(shippar.data(1,12));%
        shippar.moored.Tn3        = cell2mat(shippar.data(1,13));%
        shippar.moored.Tn5        = cell2mat(shippar.data(1,14));%
    end
        
end
%%%%%%%%%%%%%%%%%%OptionS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.interior.check   =GUImain.InteriorProp;%1 -> calculated, 0-> No   %%Interior flows properties
if options.interior.check==1
    options.interior.time    =[GUImain.IP_time(1);GUImain.IP_time(2);GUImain.IP_dt] ; %par.IP_time(1)=tstart ; par.IP_time(2)=tend; par.IP_time(3)=#dt;
end
if GUImain.partitionId==0
    options.partition.N=GUImain.partition;
    options.partition.auto='off';
    options.partition.combine=GUImain.partitionCombines;
else
    options.partition.auto='on';
    options.partition.combine=GUImain.partitionCombines;
end
options.debug        = 'No' ;  %debug the process?
options.log          = 'Yes';  % save or not log after finishing
options.OnlyEta      = 'No';  %save only surface elevation to reduce size of the MAT-file to be saved
options.IdFlagGui    = GUImain.FlagGui;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model.dyn = DynModelName(model,bath);

%% DEFAULT PARAMETERS
par.g               =  9.81;            % gravitational acceleration [m/s^2]
par.rho             =  1;               % mass density of water;
par.dispinterpol    =  2;               % # of interpolation curves for dispersion; default = only choice???
% par for ode:
par.odetol          =  1e-4;            % relative error tolerance in ODE-solver
par.odesolver       = 'ode45';          % ode-solver to be used for time integration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


PreProcessing;


