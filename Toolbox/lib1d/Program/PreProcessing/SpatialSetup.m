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
%%%%%%%%%%%%%%%%%  Spatial Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx=round((spatial.Xright-spatial.Xleft)*10000/(spatial.dx*10000))+1;
par.Nx  = Nx; 
x           = linspace(spatial.Xleft,spatial.Xright,par.Nx);
par.x       = x;
par.dx      = x(2)-x(1);
par.Xleft   = x(closest(x,spatial.Xleft));%
par.Xright  = x(closest(x,spatial.Xright));%
par.Xinflux = x(closest(x,spatial.Xinflux));
%%%%
par.bf0     = bath.FBL;
%par.dampchar= dampzone(x,par.bf0);
par.cfSA    = cf(x,par.Xleft,par.bf0(1))-cf(x,par.Xright-par.bf0(2),par.bf0(2)); % cfSimulationArea
 % cfSimulationArea

%% Fourier space (physical wavenumbers) %%%%%%%%%%%%%%%%%%%%%%%
par.k           = freqspace(x);

par.dk          = par.k(2)-par.k(1);
%% Bathymetry
if  strcmp(bath.name, 'Flat')
    par.depth   =  bath.depth;
    par.bathy   = -par.depth*ones(length(x),1);
elseif strcmp(bath.name,'Slope')
    par.bathy   = Bathy_slope(x,bath.par)';
    indXinf     = closest(x,par.Xinflux);
    par.depth   = -par.bathy(indXinf);
elseif strcmp (bath.name,'Slope (Shore)')   % for slope or run-up
    [depthX,bath.par] = Bathy_runup(x,par.Xinflux,bath.par,par.cfSA);
    par.bathy         = -depthX;
    indXinf           = closest(x,par.Xinflux);
    par.depth         = depthX(indXinf);
else % fun.bath
    par.bathy   = bathy_data(bath,x,par.cfSA);
    indXinf     = closest(x,par.Xinflux);
    par.depth   = -par.bathy(indXinf);
end
if max(any(par.bathy>=0))==1
   par.bathyplus=zeros(size(par.bathy));
   par.bathyplus(par.bathy>=0)=par.bathy(par.bathy>=0);
   par.bathymin=zeros(size(par.bathy));
   par.bathymin(par.bathy<0)=par.bathy(par.bathy<0);
   if par.bathy(1)<par.bathy(end)
   par.IndS=find(par.bathy>0,1,'first');
   else
   par.IndS=find(par.bathy>0,1,'last');
   end
       
end
%%%%%%%%%%%%%%% friction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.friction =bath.friction;
if par.friction.check==1
[par.friction.index,par.friction.Cf]=Index_friction(x,par.friction.data);
end

%%%%%%%%%%%%%%% Wall %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.wall.presence=bath.wall.check;
if strcmp(par.wall.presence,'Yes')
    par.wall.length     =bath.wall.length;
    par.wall.refl_Coef  =bath.wall.Coef;
    par.wall.alphaR     =(1-par.wall.refl_Coef)/(1+par.wall.refl_Coef);
    Xwall               =x(closest(x,bath.wall.position));
    par.wall.method    =bath.wall.method;
    par.wall.type       =bath.wall.type; 
    par.wall.position   =Xwall; 
    par.wall.Evmodes    =bath.wall.Evmodes;
    dxWallchar=zeros(size(x));ddxWallchar=zeros(size(x));
    if par.wall.type==1
        if par.wall.position > par.Xinflux
            Wallchar                =Heaviside(Xwall-x);
            indW0                   =find(Wallchar==0,1,'first');
            CutOff                  = 1-ReflCoef_to_cutoff1(par.wall.refl_Coef,model.nonlinear);
            %CutOff                 =ReflCoef_to_cutoff2(par.wall.refl_Coef); %Analytic
            Wallchar(indW0:end)     =CutOff;
            SignProp                =1;
        else
            Wallchar              = Heaviside(x-Xwall);
            indW0                 = find(Wallchar==0,1,'last');
            CutOff                = 1-ReflCoef_to_cutoff1(par.wall.refl_Coef,model.nonlinear);
            %CutOff               = ReflCoef_to_cutoff2(par.wall.refl_Coef); %Analytic
            Wallchar(1:indW0)     = CutOff;
            SignProp              = -1;
        end
        par.wall.Rho=Wallchar;
        
    else
        if  par.wall.position > par.Xinflux
            Wallchar      = Heaviside(Xwall-x);
            indW0         =find(Wallchar==0,1,'first');
            SignProp      =1;
        else
            Wallchar = Heaviside(x-Xwall);
            indW0         =find(Wallchar==0,1,'last');
            SignProp      =-1;
        end
            par.wall.char =Wallchar;
          %  par.wall.file_def= bath.wall.file_R_def;
            par.wall.file    =bath.wall.file_R;
            %%%Setup par.wall.Rho  in GenerationSetup.m because required information of frequency domain       
    end
%     if  par.wall.position>par.Xinflux
%         par.wall.wallInfdir='Uni-';
%     else
%         par.wall.wallInfdir='Uni+';
%     end
    par.wall.wallInfdir='Bi';
     
    par.wall.indW0=indW0;
    par.wall.SignProp=SignProp;
    if par.wall.method==1
        par.wall.gamX_Check='Area';
        [par.wall.gamX_hat,par.wall.gamX_skew_hat]=funWall_FluxSpatGam(par);
    end
end