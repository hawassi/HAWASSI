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

cutfrac = par.cutfrac; k= par.k; x = par.x; k_p  = influx.k_p; dtout = par.dt;
depth_inf  = influx.depth;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Om         = str2func(model.dispersion);
omAdd      = model.OmFun;
Up         = str2func(['Up',model.dispersion(3:end)]);

Oprt.Omd   = Om(k,depth_inf,omAdd);
Oprt.Upd   = Up(k,depth_inf,omAdd);
Oprt.Csq   = Oprt.Upd.^2;
Oprt.F0    = Oprt.Csq.*1i.*k/g;
Oprt.L0    = -1i*k.*Oprt.F0;
Oprt.M0    = -1i*k.*Oprt.Csq/g;
Oprt.M1    = -Oprt.Csq.*1i.*k/g; % M0=M1 in Homogeneous case
% it will be not same in
% nonHomogeneous case such as a Wall

if  input.type~=1
    Oprt.Cpeak = Om(k_p,depth_inf,omAdd)/k_p;
    CpLeft=Om(k_p,-par.bathy(1),omAdd)/k_p;
    CpRight= Om(k_p,-par.bathy(end),omAdd)/k_p;
else
    Oprt.Cpeak =sqrt(g*depth_inf(1));
    CpLeft= sqrt(g*-par.bathy(1));
    CpRight= sqrt(g*-par.bathy(end));
end
Oprt.aal         = fun_alias(k,cutfrac);

par.dampcoef=ones(length(par.x),1);
par.dampcoef(1:floor(end/2),1)=7*CpLeft/(0.9*par.bf0(1)); % factor 0.9 to damp faster before end of the domain
par.dampcoef(floor(end/2)+1:end,1)=7*CpRight/(0.9*par.bf0(2));

if shippar.check==1 && (strcmp(bath.type,'F')||strcmp(bath.type,'B') || strcmp(bath.type,'U'))
    Sshape=shippar.form.shapeXcZc(end,:).';
    DD=-par.bathy+Sshape;
    D_mid= (min(DD)+max(DD))/2;
    
    if strcmp(bath.type,'B') || strcmp(bath.type,'U')
        Dmin=min(-par.bathy)+min(shippar.form.shapeXcZc(end,:).');
        D_mid=mean(DD(Sshape<0));
    else
        Dmin=min(DD);
    end
    Dplus=max(DD);
   
    DD=linspace(Dmin,Dplus,100);
    [D_min,D_plus,~,~,~,Spgammin,Spgammid,Spgamplus] = Up3IP(k,Up,Om,omAdd,influx.nu_p,DD,D_mid,...
        Proj.Dir,Proj.savename,model.dyn);
    Oprt.S.Csqmin      = Up(k,D_min,omAdd).^2;
    Oprt.S.Csqmid      = Up(k,D_mid,omAdd).^2;
    Oprt.S.Csqplus     = Up(k,D_plus,omAdd).^2;
    Oprt.S.Sp_gammin      = Spgammin;
    Oprt.S.Sp_gammid      = Spgammid;
    Oprt.S.Sp_gamplus     = Spgamplus;
    Oprt.S.interp      =3;
end

if strcmp(bath.type,'B') || strcmp(bath.type,'U')
    DD=-par.bathy;
    if bath.interpRef==2
        [D_min,D_plus,gammin,gamplus] = Up2IP(k,Up,Om,omAdd,influx.nu_p,DD,...
            Proj.Dir,Proj.savename,model.dyn);
        Oprt.Csqmin      = Up(k,D_min,omAdd).^2;
        Oprt.Csqplus     = Up(k,D_plus,omAdd).^2;
        Oprt.gammin      = gammin;
        Oprt.gamplus     = gamplus;
        Oprt.interp      = 2;
        
        Oprt.Ommin       = Om(k,D_min,omAdd);
        Oprt.Omplus      = Om(k,D_plus,omAdd);
        
    else
        D_mid= bath.DmidRef;
        [D_min,D_plus,gammin,gammid,gamplus,Spgammin,Spgammid,Spgamplus] = Up3IP(k,Up,Om,omAdd,influx.nu_p,DD,D_mid,...
            Proj.Dir,Proj.savename,model.dyn);
        Oprt.Ommin       = Om(k,D_min,omAdd);
        Oprt.Ommid       = Om(k,D_mid,omAdd);
        Oprt.Omplus      = Om(k,D_plus,omAdd);
        
        Oprt.Csqmin      = Up(k,D_min,omAdd).^2;
        Oprt.Csqmid      = Up(k,D_mid,omAdd).^2;
        Oprt.Csqplus     = Up(k,D_plus,omAdd).^2;
        Oprt.gammin      = gammin;
        Oprt.gammid      = gammid;
        Oprt.gamplus     = gamplus;
        Oprt.interp      =3;
        
        %     Oprt.Sp_gammin      = Spgammin;
        %     Oprt.Sp_gammid      = Spgammid;
        %     Oprt.Sp_gamplus     = Spgamplus;
    end
    
    
    
elseif strcmp(bath.type,'BR') || strcmp(bath.type,'UR')
    % % if input.type~=1
    % % maxEtaInit  = max(influx.timesig_insig(:,2));
    % % else
    % % maxEtaInit =max(Ifft(IVP.zeta0_hat));
    % % end
    maxEtaInit=1.5.*influx.Hs;
    detaC2=str2func(model.detaC2); %for Hsdirect model
    Cder=str2func(model.Cder);
    
    H_mid=bath.DmidRef;
%     if shippar.check==1
%         Depth=-par.bathy+shippar.form.shapeXcZc(end,:).';
%     else
        Depth=-par.bathy;
%    end
%     if model.nonlinear==2
%         H_max = max(Depth(1),Depth(end))+maxEtaInit;
%     elseif model.nonlinear==3
        H_max = max(Depth(1),Depth(end));
%    end
    k_cut=max(k)/par.cutfrac;
    H_minShore=bath.hmin; %(3*influx.nu_p/4/k_cut)^2/par.g/2;
    H_min=5*H_minShore;
    %H_min  = min(H_plus/100,1);%0.8;% %~~0.1 for liu case%
    
    %%%%%%%%%%%%%%%%%%%%%rhsshore
    if model.nonlinear>=2 %200917 Now second order also use polynomial model
        Dperlambda=1/200; %% this is for very longwaves i.e T=100 s, if 1/20 (SWE) this leads to D~~200 (too much)
        KD=2*pi*Dperlambda;
        H_min=max(H_min,par.g*(KD*tanh(KD))/(influx.nu_p.^2)); % for windwave H_min=5*Hminshore it seems safe.
    end
    
    [IntCoef,H_min1,H_max1,H_mid1,H_min2,H_max2,H_mid2]= Up3IP_2part_runup(k,Up,Om,omAdd,influx.nu_p,par.cutfrac,g,Depth,H_min,H_mid,...
        H_max,Proj.Dir,Proj.savename,model,influx.categoryMaxDepth);
    % [IntCoefdetaC2]= Up2Derv3IP_runup(k,Up,Om,omAdd,detaC2,influx.nu_p,g,Depth,H_min,H_mid,...
    %         H_max,Proj.Dir,Proj.savename,model);
    
    % [IntCoefdetaC2]= detaC23IP_2part_runup(k,Up,Om,omAdd,detaC2,influx.nu_p,g,Depth,H_min,H_mid,...
    %         H_max,Proj.Dir,Proj.savename,model);
    par.interp.IntCoef=IntCoef;
    %par.interp.IntCoefdetaC2=IntCoefdetaC2;
    
    par.interp.H_min1  =H_min1;
    par.interp.H_mid1  =H_mid1;
    par.interp.H_max1  =H_max1;
    par.interp.H_min2  =H_min2;
    par.interp.H_mid2  =H_mid2;
    par.interp.H_max2  =H_max2;
    
    par.interp.H_min  =H_min;
    par.interp.H_mid  =H_mid;
    par.interp.H_max  =H_max;
    
    
    par.interp.H_minShore=H_minShore;
    
    Oprt.C2m1           =Up(k,H_min1,omAdd).^2;
    Oprt.C2p1           =Up(k,H_max1,omAdd).^2;
    Oprt.C2c1           =Up(k,H_mid1,omAdd).^2;
    Oprt.C2m2           =Up(k,H_min2,omAdd).^2;
    Oprt.C2p2           =Up(k,H_max2,omAdd).^2;
    Oprt.C2c2           =Up(k,H_mid2,omAdd).^2;
    
    if  input.type~=1  %for damping coefficient
        Oprt.CpeakL = Om(k_p,Depth(1),omAdd)/k_p;
        Oprt.CpeakR = Om(k_p,Depth(end),omAdd)/k_p;
    else
        Oprt.CpeakL =sqrt(g*Depth(1));
        Oprt.CpeakR =sqrt(g*Depth(end));
    end
    
    par.dampcoef=ones(length(par.x),1);
end

