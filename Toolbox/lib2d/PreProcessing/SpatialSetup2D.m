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
%%%%%%%%%% Spatial Setup                                         %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Physical Space
Nx=round((spatial.dom.x(2)-spatial.dom.x(1))*10000/(spatial.dom.dx*10000))+1;
dom.Nx  = Nx; 
dom.X   = linspace(spatial.dom.x(1),spatial.dom.x(2),dom.Nx);
dom.dx  = dom.X(2)-dom.X(1);
Ny=round((spatial.dom.y(2)-spatial.dom.y(1))*10000/(spatial.dom.dy*10000))+1;
dom.Ny  = Ny; 
dom.Y   = linspace(spatial.dom.y(1),spatial.dom.y(2),dom.Ny);
dom.dy  = dom.Y(2)-dom.Y(1);
[dom.XX,dom.YY]  = meshgrid(dom.X,dom.Y);

%% Fourier space (physical wavenumbers) %%%%%%%%%%%%%%%%%%%%%%%
dom.kx           = funC_freqspace(dom.X)';
dom.dkx          = dom.kx(2)-dom.kx(1);
dom.ky           = funC_freqspace(dom.Y)';
dom.dky          = dom.ky(2)-dom.ky(1);
[dom.Kx,dom.Ky]  = meshgrid(dom.kx,dom.ky);
dom.KK           = sqrt(dom.Kx.^2+dom.Ky.^2);
dom.cutfracx     = spatial.dom.cutfracx;
dom.cutfracy     = spatial.dom.cutfracy;
dom.aal          = funSP_alias2d(dom); %wavenumbercut (anti-alias)
%assignin('base','dom',dom);
%%%%%%%%%%%%%%%%% FBdy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(spatial.fbl.option,'Yes')
[dom.cfSA,dom.fbl]=funD_simulationAreaChar(dom,spatial.fbl);%funC_cfSA2d(dom.X,dom.Y,dom.fbl);  
else
dom.fbl.l=0;dom.fbl.r=0;
dom.fbl.b=0;dom.fbl.t=0;
dom.cfSA=ones(size(dom.XX));
dom.fbl.char=zeros(size(dom.XX));
end
%%%%%%%%%%%%%%  Wall   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dom.wall=input.wall;
if strcmpi(dom.wall.option,'yes')
  % [dom.wall.char,dom.wall.gradchar]= funW_characteristic(dom,model.nonlinear);   
    [dom.wall.charAll,dom.wall.char,dom.wall.gradchar,dom.wall.charInfl,dom.wall.NInfl,dom.wall.ReflCoef]= ...
       funW_characteristic_for_combined_method(dom,model.nonlinear);%funW_characteristic_for_infl_method(dom);   
   
   dom.wall.bdyInfl.index   = funW_indx_bdyCharInfl(dom.wall.charAll,dom.wall.charInfl);  
   dom.wall.bdyCutK.index1   = funC_indx1_bdy(dom.wall.char);  
   dom.wall.bdyCutK.index0   = funC_indx0_bdy(dom.wall.char);  
   
    dom.wall.N          = length(dom.wall.param(:,1));
%    FBLdat=[dom.XX(dom.wall.bdy.index) dom.YY(dom.wall.bdy.index)]; 
%    save FBLdat FBLdat
   dom.wall.bdyInfl.dirflag = funW_bdy_walldir_flag(dom);
        
   if any(dom.wall.char<=0.05)
        dom.fbl.char(dom.wall.char<=0.05)=1./(0.9*sqrt(dom.dx^2+dom.dy^2));  
        dom.cfSA(dom.wall.char<=0.05)=0;
%         figure;
%         subplot(2,1,1)
%         surf(dom.XX,dom.YY,dom.wall.char,'edgecolor','none');view(2);
%         subplot(2,1,2)
%         surf(dom.XX,dom.YY,dom.fbl.char,'edgecolor','none');view(2);        
   end
%    MMM=zeros(size(dom.XX));
%    MMM(dom.wall.bdyCutK.index0)=1;
%    [X,Y]=meshgrid(dom.X,dom.Y);
%    ff=figure(104);
%    set(ff,'Renderer','zbuffer'); %due to graphics driver
%    surf(X,Y,MMM,'Edgecolor','none');
%    view(2);
else
    dom.wall.char = 1; %% added by Nida
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bathymetry
if  strcmpi(bath.name, 'Flat')
    dom.bathy.profile   = -bath.depth*ones(dom.Ny,dom.Nx);
elseif strcmpi(bath.name,'Slope')
    dom.bathy.profile   = funB_Bathy_slope(dom,bath);
    dom.bathy.interp    = bath.interp;
elseif strcmpi(bath.name,'Shore')
    dom.bathy.profile   = funB_Bathy_shore(dom,bath);
    dom.bathy.interp    =bath.interp;
elseif strcmpi(bath.name,'User-defined')
    dom.bathy.profile   = funB_Bathy_userdefined(dom,bath);
    dom.bathy.interp    = bath.interp;
    ID_shore=0;
    if any(dom.bathy.profile>0,1) 
        ID_shore=1;
    end
    if any(dom.bathy.profile>0,2) 
        ID_shore=1;
    end
    
    if ID_shore==1
        disp('here')
        bath.type=funB_userdefined_shoretype(dom.bathy.profile,dom.XX,dom.YY,bath.Hmin);
        bath.name='Shore';
        dom.bathy.profile=funB_smoothed_bathy_userdefined_runup(dom,bath.type);
    end
end

if strcmpi(bath.name,'Shore')
   dom.bathy.plus=dom.bathy.profile;dom.bathy.plus(dom.bathy.plus<0)=0;
   dom.bathy.min=dom.bathy.profile;dom.bathy.min(dom.bathy.min>0)=0; 
end
%%%%%%%%%%%%%%% friction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dom.friction.input =bath.friction;
if dom.friction.input.check==1
[dom.friction.Char]=funB_frictionChar2D(dom,bath.friction);
end

%%%%%%%%%%%%%%% Floating body %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
body=input.body;
if body.option     == 1
[body.par]= parambodysetup(dom,input.body);
end