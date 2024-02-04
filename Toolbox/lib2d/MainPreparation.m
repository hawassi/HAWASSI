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
%%%%%%%%%% version: 24 October 2016                                  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RK%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Main-User                                             %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------------------------------------------------------%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
%------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%Project Preparation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Preproc=MainPreparation(GUIinput,jProgressBar,statusbarObj)
Proj.path                = GUIinput.proj.path;
Proj.savename            = GUIinput.proj.name;
Proj.usernote            = GUIinput.proj.note;
Proj.projdir             = GUIinput.proj.projdir;
Proj.workdir             = GUIinput.proj.workdir;%
Proj.module              = GUIinput.proj.module;%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model.evol             =  GUIinput.modeldyn(1:2);  % 'HS' for dispersive Ham. System
model.nonlinear        =  str2num(GUIinput.modeldyn(3));%GUImain.nonlin;      % for nonlin order: #=1(linear),2,3,4;
model.breaking.check   =  GUIinput.modelbreak.check;%GUImain.Break;       % 'Yes' or 'No'
if strcmpi(model.breaking.check,'Yes')
    model.breaking.KBC  = GUIinput.modelbreak.initiation;    % Kin.Br.Cond.: U/C in [0.7:1]
    model.breaking.TC   = GUIinput.modelbreak.termination;
    model.breaking.Tchar= GUIinput.modelbreak.Tstar;
    model.breaking.delb    = 1.5;              % Mixing length Coef
end
model.dispersion       =  ['Om2d',GUIinput.modeldisp];%GUImain.disp;        % OmExact,OmSWE,OmKdV,OmBBM,  available; specify OmUser
model.phasevel         =  ['Cp2d',GUIinput.modeldisp];
model.groupvel         =  ['Ug2d',GUIinput.modeldisp];%GUImain.groupvel;    % specify UgUser related to OmUser
model.OmFun            =  '';%GUImain.dispEqInput;
model.UgFun            =  '';%GUImain.groupVelEqInput;

model.current.check   =GUIinput.modelcurrent.check;
model.current.ux      =GUIinput.modelcurrent.ux;
model.current.uy      =GUIinput.modelcurrent.uy;
model.phiForm         =GUIinput.modelphiForm  ;% 0 in u,v formulation, 
                                                % 1 in phi formulation.
                                      
model.runupevo        ='Hs'; % Hs or Hsdirect;
if model.phiForm      ==1
model.runupevo        ='Hsdirect'; % Hs or Hsdirect;    
end
                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%Initial Conditions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ivp.typename           = GUIinput.wave_ivp.type;        %GUImain.IVP.typename;%#1: Zero, 2:Gaussian, 3: N-Wave, 4:User-defined 
ivp.A                  = GUIinput.wave_ivp.A;       %GUImain.IVP.A;
ivp.sigma              = GUIinput.wave_ivp.stdev;       %GUImain.IVP.sigma; %standard deviation
ivp.x0                 = GUIinput.wave_ivp.centerposition;       %GUImain.IVP.x0;
ivp.userdata           = GUIinput.wave_ivp.userdata;       %GUImain.IVP.file.data;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%Wave-Input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input.wave.option          = GUIinput.wave_influx.type  ; 
input.wave.propertiesdata  = GUIinput.wave_influx.propertiesdata;
input.wave.methoddata      = GUIinput.wave_influx.methoddata;
input.wave.userdata        = GUIinput.wave_influx.userdata;

input.wave.ramp.check      = GUIinput.wave_influx.rampcheck;            %GUImain.ramp.check;   %
input.wave.ramp.lengthfact = GUIinput.wave_influx.rampfactor;            %GUImain.ramp.val;     % the length of the ramp function #Tp
input.wave.adjzone.check   = GUIinput.wave_influx.nonlinadjcheck; %#*lambda_p
input.wave.adjzone.lengthfact=GUIinput.wave_influx.nonlinadjfactor;
input.wave.tapered.check   =GUIinput.wave_influx.ramplinecheck;
input.wave.tapered.length  =GUIinput.wave_influx.ramplinefactor/2;%   %each side of line taperd by a half of the fraction of the length of the line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% Boudnary assimilation Wave-input%%%%%%%%%%%%%%%%%%%% 
if GUIinput.wave_bdy_assim.checkVal==1
input.bdyassim.option=0;
input.bdyassim.check='No';
else
input.bdyassim.option=1;   
input.bdyassim.check='Yes';
end
input.bdyassim.shapeOpt = GUIinput.wave_bdy_assim.shapeCheck;
input.bdyassim.shapedata=GUIinput.wave_bdy_assim.shapeUserdata;
input.bdyassim.halfcirc_R1 = GUIinput.wave_bdy_assim.R1;
input.bdyassim.halfcirc_xc = GUIinput.wave_bdy_assim.xc;
input.bdyassim.halfcirc_yc = GUIinput.wave_bdy_assim.yc;
input.bdyassim.smoothfact=GUIinput.wave_bdy_assim.smoothfact;
input.bdyassim.assimdata=GUIinput.wave_bdy_assim.assimdata;
if model.phiForm==1
input.bdyassim.assimdata_phi=GUIinput.wave_bdy_assim.assimdata_phi;
input.bdyassim.cb_phi=GUIinput.wave_bdy_assim.cb_phi;
else
input.bdyassim.assimdata_u=GUIinput.wave_bdy_assim.assimdata_u;
input.bdyassim.assimdata_v=GUIinput.wave_bdy_assim.assimdata_v;
input.bdyassim.cb_vel=GUIinput.wave_bdy_assim.cb_vel;    
end
input.bdyassim.cb_nonlinAdj=GUIinput.wave_bdy_assim.cb_nonlinAdj;
input.bdyassim.nonlinAdj_distance=GUIinput.wave_bdy_assim.nonlinAdj_distance;
input.bdyassim.nonlinAdj_smooth=GUIinput.wave_bdy_assim.nonlinAdj_smooth;

input.bdyassim.propdir=GUIinput.wave_bdy_assim.propdir;
input.bdyassim.tinit=GUIinput.wave_bdy_assim.tinit;
input.bdyassim.tend=GUIinput.wave_bdy_assim.tend;
input.bdyassim.dt=GUIinput.wave_bdy_assim.dt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Time simulation setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeSimul.t_init=GUIinput.time.t_start;
timeSimul.t_end=GUIinput.time.t_end;
timeSimul.dt=GUIinput.time.dt;
%%%%%%%%%  Spatial domain setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
spatial.dom.x                 =[GUIinput.spatialinterv.xmin GUIinput.spatialinterv.xmax]; %[min max]
spatial.dom.y                 =[GUIinput.spatialinterv.ymin GUIinput.spatialinterv.ymax]; %[min max]
spatial.dom.dx                =GUIinput.spatialgrid.dx; %[min max]
spatial.dom.dy                =GUIinput.spatialgrid.dy; %[min max]
spatial.dom.cutfracx          =GUIinput.fourier_cutfrac.k;
spatial.dom.cutfracy          =GUIinput.fourier_cutfrac.k;

%Fourier Boundary 
spatial.fbl.option  =GUIinput.spatialdamp.check;
spatial.fbl.param   =GUIinput.spatialdamp.data;
spatial.fbl.userdata=GUIinput.spatialdamp.userdata;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%---Wall---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input.wall.option        =GUIinput.spatialwall.check;
input.wall.param         =GUIinput.spatialwall.data;  
input.wall.userdata      =GUIinput.spatialwall.userdata;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%Ship/Struct-Input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input.body.option        = 0;%             %#1=Include ship, #0=Not include 
% input.body.fixed         = 1;%             %#1=Fixed, #0=Free
% if input.body.option     == 1;
% input.body.N             = 1;      %number of structure
% input.body.RBdyn         ='heave';
% input.body.form          ={'vertical cylinder'};
% if strcmp(input.body.form,'barge')
% input.body.length        =0.5;
% input.body.draft         =0.25;
% input.body.width         =2.76;  %
% input.body.xyz0         =[10 1.38 0];
% elseif strcmp(input.body.form,'vertical cylinder')
% input.body.draft         =0.6;
% input.body.radius        =0.15;
% input.body.xyz0         =[4.375 0.6 0];
% end
% input.body.interpOprt    =3; %#2 or 3 points interpolation
% input.body.saveoutputvar =1; %#Save ode var 1, not 0
% end

%%%%%%%%%%%%%%%%%%%%%%%%Bathymetry%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bath.type               =GUIinput.spatialbathy.type;

if strcmpi(bath.type,'Flat')==1
    bath.depth          = GUIinput.spatialbathy.depth; %GUImain.depth;       %depth in [m]
    bath.name           = 'Flat';
elseif strcmpi(bath.type,'Slope (in x-axis)')||strcmpi(bath.type,'Slope (in y-axis)')
    bath.par            = [GUIinput.spatialbathy.xymindepth GUIinput.spatialbathy.xymaxdepth...
        GUIinput.spatialbathy.slope GUIinput.spatialbathy.startslope];%[GUImain.depth(1);GUImain.depth(2);GUImain.slope(1);GUImain.slope(2)];%[30,10,1/10,200]; % for slope: par.bath = [Depth_left,Depth_right,slope,X_startslope];
    bath.name           = 'Slope';
    bath.interp         = GUIinput.spatialbathy.interp;
    if GUIinput.spatialbathy.interp > 2 % nunu
        bath.depthref   = GUIinput.spatialbathy.xymiddepth;
    else
        bath.depthref   = [];
    end
elseif strcmpi(bath.type,'Shore (in x-axis)')||strcmpi(bath.type,'Shore (in y-axis)')          %bathy runup
    bath.par            = [GUIinput.spatialbathy.maxdepthshore ...
    GUIinput.spatialbathy.slopeshore GUIinput.spatialbathy.shoreposition];%[GUImain.depth;GUImain.slope(1);GUImain.slope(2)];
    bath.name           = 'Shore';
    bath.interp         = GUIinput.spatialbathy.interp;
    if GUIinput.spatialbathy.interp > 2 % nunu
        bath.depthref   = GUIinput.spatialbathy.xymiddepth;
    else
        bath.depthref   = [];
    end
    bath.Hmin           =GUIinput.spatialbathy.mindepthshore;
else
    %load('D:\Workspace\1. Codes\developer\HaWaSSI\HAWASSI_AB2_160809_Berkhof_Exp\Testcases\Berkhoff_82\BathyBerkhoff82.mat')
    bath.name           = 'User-defined';
    bath.userdata       = GUIinput.spatialbathy.userdata.data;%GUImain.bathyO.data;% 'Bathymetry data'
    bath.interp         = GUIinput.spatialbathy.interp;
    if GUIinput.spatialbathy.interp > 2 % nunu
        bath.depthref   = GUIinput.spatialbathy.xymiddepth;
    else
        bath.depthref   = [];
    end
   
    if any(bath.userdata(:,3)>0)
       bath.Hmin           =GUIinput.spatialbathy.mindepthshore;
    end
end
bath.friction.check       = GUIinput.spatialfriction.check;
bath.friction.param      = GUIinput.spatialfriction.data;
bath.friction.userdata   = GUIinput.spatialfriction.userdata;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%OptionS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.interior.check   =GUIinput.option_intflow.check; 
if options.interior.check==1
options.interior.time    =[GUIinput.option_intflow.tinit GUIinput.option_intflow.tend ...
    GUIinput.option_intflow.dt_fact];
end
options.partitionCheck   =GUIinput.option_partition.check_def;
options.partitionCombines=GUIinput.option_partition.check_combine;
if options.partitionCheck==0
options.partition.N=GUIinput.option_partition.total;
options.partition.auto='off';
options.partition.combine=options.partitionCombines;
else
options.partition.auto='on'; 
options.partition.combine=options.partitionCombines;
end
options.mc.check=GUIinput.option_mc.check;
options.mc.waveinput_check=GUIinput.option_mc.waveinput_check;
options.mc.waveinput_data=GUIinput.option_mc.waveinput_userdata;
options.mc.numrun_check=GUIinput.option_mc.numrun_check;
options.mc.numrun=GUIinput.option_mc_numrun_edit;

options.mc.numset_check=GUIinput.option_mc.numset_check;
options.mc.numset_px=GUIinput.option_mc.numset_px;
options.mc.numset_py=GUIinput.option_mc.numset_py;
options.mc.N=funP_N_mc(options.mc);

options.debug        = 'No' ;  %debug the process? 
options.log          = 'Yes';  % save or not log after finishing 
if GUIinput.option_outputvar==1
options.OnlyEta      = 'Yes';   %save only surface elevation to reduce size of the MAT-file to be saved
else
options.OnlyEta      = 'No';   %save only surface elevation to reduce size of the MAT-file to be saved
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model.dyn = funP_DynModelName(model,bath);    

%% DEFAULT PARAMETERS
par.g               =  9.81;            % gravitational acceleration [m/s^2]
par.rho             =  1;               % mass density of water;
% par for ode:
par.ode.tol          =  GUIinput.option_ode_tol; % relative error tolerance in ODE-solver

if GUIinput.option_ode_sol==1
par.ode.solver       = 'ode45';          % ode-solver to be used for time integration 
elseif GUIinput.option_ode_sol==2
par.ode.solver       = 'ode23';    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PreProcessing2D;
end
% if Execute==2||Execute==12
% ODEsolver;
% end
