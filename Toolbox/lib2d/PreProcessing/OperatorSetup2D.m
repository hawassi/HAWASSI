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
g            =par.g;
if model.current.check==0
Oprt_Om2d_fun= str2func(['funOprt_',model.dispersion]);
omAdd        = model.OmFun;
Oprt_Up2d_fun= str2func(['funOprt_',model.phasevel]);
else
Oprt_Om2d_fun= str2func(['funOprt_',model.dispersion,'_current']);
omAdd        = model.OmFun;
Oprt_Up2d_fun= str2func(['funOprt_',model.phasevel,'_current']);    
end


Om_fun       = str2func(['funOprt_Om',model.dispersion(5:end)]); 
Up_fun       = str2func(['funOprt_Up',model.dispersion(5:end)]);%@(k,d,g,omAdd) Om_fun(k,d,g,omAdd)./k;
Ug_fun       = str2func(['funOprt_Ug',model.dispersion(5:end)]);

% nunu
Depth_all    = -dom.bathy.profile(:);
Oprt.Cpeak   = Up_fun(min(k_p),meandepth,g,omAdd);
Oprt.Cpeak_A = reshape(Up_fun(min(k_p),Depth_all,g,omAdd),dom.Ny,dom.Nx);	

if strcmp(bath.type,'Flat') && input.body.option==0
    if model.current.check==0
    Oprt.Om2d  =Oprt_Om2d_fun(dom.Kx,dom.Ky,meandepth,g,[1 1]);
    else
    Oprt.Om2d  =Oprt_Om2d_fun(dom.Kx,dom.Ky,meandepth,g,[1 1],model.current.ux,...
        model.current.uy);    
    end
    Oprt.L2d   =Oprt.Om2d.^2/g;
     if model.current.check==0
        Oprt.Cp2d =Oprt_Up2d_fun(dom.Kx,dom.Ky,meandepth,g);
    else
        Oprt.Cp2d =Oprt_Up2d_fun(dom.Kx,dom.Ky,meandepth,g,model.current.ux,...
            model.current.uy);
    end
    Oprt.Cp2dSq=Oprt.Cp2d.^2;
    
elseif ~strcmpi(bath.type,'Flat') && input.body.option==0
    if strcmpi(bath.name,'Shore') %% Nida
        [InterpDisp,dom.bathy.HminShore,dom.bathy.HminDisp]= funOprt_DispInterpolation_3p_runup(g,influx,dom,Ug_fun,Up_fun,Om_fun,omAdd,nupeak,...
            -dom.bathy.profile,bath,Proj,model.runupevo);
        Oprt.Cp2dSqInterp = funOprt_Cp2dSq_InterpComponents_3p_runup(g,bath,...
            Oprt_Up2d_fun,InterpDisp,dom,model.current);
       if strcmpi(model.runupevo,'Hs')
        Oprt.Cp2dInterp = funOprt_Cp2d_InterpComponents_3p_runup(g,bath,...
            Oprt_Up2d_fun,InterpDisp,dom,model.current);
        Oprt.Om2dSqInterp = funOprt_Om2dSq_InterpComponents_3p_runup(g,bath,...
            Oprt_Om2d_fun,InterpDisp,dom,model.current);
       end
   else
        if dom.bathy.interp==3
            [InterpDisp]= funOprt_DispInterpolation_3p(g,dom,Ug_fun,Up_fun,Om_fun,omAdd,nupeak,...
                -dom.bathy.profile,bath,Proj,model.dyn);
            Oprt.Om2dSqInterp=funOprt_Om2dSq_InterpComponents_3p(g,Oprt_Om2d_fun,InterpDisp,dom,model.current);
            Oprt.Cp2dSqInterp =funOprt_Cp2dSq_InterpComponents_3p(g,Oprt_Up2d_fun,InterpDisp,dom,model.current);
            if strcmpi(dom.wall.option,'Yes')
                Oprt.Cp2dInterp =funOprt_Cp2d_InterpComponents_3p(g,Oprt_Up2d_fun,InterpDisp,dom,model.current);
               % Oprt.Cp2dSqInterp =funOprt_Cp2dSq_InterpComponents_3p(g,Oprt_Up2d_fun,InterpDisp,dom,model.current);
            end
            
        else
            [InterpDisp]= funOprt_DispInterpolation_2p(g,dom,Up_fun,Om_fun,omAdd,nupeak,...
                -dom.bathy.profile,bath,Proj,model.dyn);
            Oprt.Om2dSqInterp=funOprt_Om2dSq_InterpComponents_2p(g,Oprt_Om2d_fun,InterpDisp,dom,model.current);
            Oprt.Cp2dSqInterp =funOprt_Cp2dSq_InterpComponents_2p(g,Oprt_Up2d_fun,InterpDisp,dom,model.current);
           
            if strcmpi(dom.wall.option,'Yes')
                Oprt.Cp2dInterp =funOprt_Cp2d_InterpComponents_2p(g,Oprt_Up2d_fun,InterpDisp,dom,model.current);
              %  Oprt.Cp2dSqInterp =funOprt_Cp2dSq_InterpComponents_2p(g,Oprt_Up2d_fun,InterpDisp,dom,model.current);
            end
        end
    end
end

if body.option== 1
  if input.body.interpOprt==3
  [InterpDisp]= funOprt_DispBodyInterpolation_3p(g,dom,Ug_fun,Up_fun,Om_fun,omAdd,nupeak,...
                            -dom.bathy.profile+body.par.shape,bath,Proj,model.dyn);  
  Oprt.Om2dSqInterp=funOprt_Om2dSq_InterpComponents_3p(g,Oprt_Om2d_fun,InterpDisp,dom,model.current);
  Oprt.OneperOm2dSqInterp=funOprt_OneperOm2dSq_InterpComponents_3p(g,Oprt_Om2d_fun,InterpDisp,dom,model.current);
   else
  [InterpDisp]= funOprt_DispBodyInterpolation_2p(g,dom,Up_fun,Om_fun,omAdd,nupeak,...
                            -dom.bathy.profile+body.par.shape,bath,Proj,model.dyn);
  Oprt.Om2dSqInterp=funOprt_Om2dSq_InterpComponents_2p(g,Oprt_Om2d_fun,InterpDisp,dom,model.current);
  Oprt.OneperOm2dSqInterp=funOprt_OneperOm2dSq_InterpComponents_2p(g,Oprt_Om2d_fun,InterpDisp,dom,model.current);
  end
  if model.current.check==0
  Oprt.Cp2dSq =Oprt_Up2d_fun(dom.Kx,dom.Ky,influx.par.meandepth,g).^2;
  Oprt.Om2d  =Oprt_Om2d_fun(dom.Kx,dom.Ky,influx.par.meandepth,g,[1 1]);
  else
  Oprt.Cp2dSq =Oprt_Up2d_fun(dom.Kx,dom.Ky,influx.par.meandepth,g,model.current.ux,model.current.uy).^2;
  Oprt.Om2d  =Oprt_Om2d_fun(dom.Kx,dom.Ky,influx.par.meandepth,g,[1 1],model.current.ux,model.current.uy);    
  end
  Oprt.L2d   =Oprt.Om2d.^2/par.g;
  
  if body.fixed==0
  [body.radpot,body.addedmass]=RadiationPotential2dbyMinimization(g,dom,Oprt,body);
  else
  body.radpot=0;body.addedmass=0;    
  end
end

