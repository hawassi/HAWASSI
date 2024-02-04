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
%%%%%%%%%% ODe output processor                                %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global dataCrestBreak dataBreak_nodes  dataSpeedBreak bodysavevar Idstop ...
    ITERbn ITERdt



if Npartition>1
    if Idstop==0
        if time(end)<tintval(end)  %time step failure
            Idstop=1;
            Flag_timestep_failure=1;
        end
    end
end

if Idstop==2
    Idstop=1;
    Flag_timestep_failure=1;
end


if  I==Npartition || Idstop==1
    set(jProgressBar,'Maximum',timeSimul.Nt, 'Value',timeSimul.Nt);
    jProgressBar.setStringPainted( true );
    statusbarObj.setText('ODE calculation done.');
    Timeodesolver =cput;
    CompRel = Timeodesolver/DeltaTimesimul;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Producing output of ode-calculation%
%spaceinterpnum =1;
% interpolate the numerical solution (for plotting purposes);
% the smallest wavelength gets 'spaceinterpnum' points
% switch dynmodel
Nend        = length(time);
if body.option==0
     if model.phiForm==1 
     nz          = length(z_hat(1,:))/2;
     else
     nz          = length(z_hat(1,:))/3;
     end
else
    Nbody       =body.N;
    nz          = (length(z_hat(1,:))-2*6*Nbody)/2;
    RB          = zeros(Nend,6*Nbody);
end
eta         = zeros(Nend,dom.Ny,dom.Nx);
if model.phiForm==1 
phi         = zeros(Nend,dom.Ny,dom.Nx);
else
u           = zeros(Nend,dom.Ny,dom.Nx);   
v           = zeros(Nend,dom.Ny,dom.Nx);   
end

for j = 1 : length(time)
    eta2hat    = z_hat(j,1:nz);
    eta2hat    = reshape(eta2hat.',dom.Ny,[]);
    eta2temp =funC_ifft2(eta2hat);
    if model.phiForm==1 
    phi2hat    = z_hat(j,nz+1:2*nz);
    phi2hat    = reshape(phi2hat.',dom.Ny,[]);
    phi2temp   = funC_ifft2(phi2hat);
    else
    u2hat    = z_hat(j,nz+1:2*nz);
    u2hat    = reshape(u2hat.',dom.Ny,[]);
    u2temp   = funC_ifft2(u2hat);   
    
    v2hat    = z_hat(j,2*nz+1:3*nz);
    v2hat    = reshape(v2hat.',dom.Ny,[]);
    v2temp   = funC_ifft2(v2hat); 
    end
    
    if strcmpi(dom.wall.option,'yes') && dom.wall.NInfl>0
        if dom.wall.ReflCoef.FlagFreqDep==1
            for JJwall=1:dom.wall.N
                if dom.wall.ReflCoef.FreqDep(JJwall).flag==1
                    Refl_Eta=dom.wall.ReflCoef.FreqDep(JJwall).char...
                        .*funC_ifft2(eta2hat.*(dom.wall.ReflCoef.fun_K(JJwall).coef+1));
                    if model.phiForm==1
                        Refl_Phi=dom.wall.ReflCoef.FreqDep(JJ).char...
                            .*funC_ifft2(phi2hat.*(dom.wall.ReflCoef.fun_K(JJwall).coef+1));
                    else
                        Refl_u=dom.wall.ReflCoef.FreqDep(JJ).char...
                            .*funC_ifft2(u2hat.*(dom.wall.ReflCoef.fun_K(JJwall).coef+1));
                        Refl_v=dom.wall.ReflCoef.FreqDep(JJ).char...
                            .*funC_ifft2(v2hat.*(dom.wall.ReflCoef.fun_K(JJwall).coef+1));
                    end
                end
            end
        else
            Refl_Eta=(1+dom.wall.ReflCoef.uniform).*eta2temp;
            if model.phiForm==1
            Refl_Phi=(1+dom.wall.ReflCoef.uniform).*phi2temp;
            else
            Refl_u=(1+dom.wall.ReflCoef.uniform).*u2temp;  
            Refl_v=(1+dom.wall.ReflCoef.uniform).*v2temp;  
            end
        end
        eta2temp(dom.wall.bdyInfl.index)=Refl_Eta(dom.wall.bdyInfl.index);
         if model.phiForm==1
         phi2temp(dom.wall.bdyInfl.index)=Refl_Phi(dom.wall.bdyInfl.index);
         else
         u2temp(dom.wall.bdyInfl.index)=Refl_u(dom.wall.bdyInfl.index);  
         v2temp(dom.wall.bdyInfl.index)=Refl_v(dom.wall.bdyInfl.index);  
         end
    end
    if (strcmpi(bath.name,'Shore')) %% added by Nida for etainit0
        eta2temp(dom.bathy.profile>0) = eta2temp(dom.bathy.profile>0)+dom.bathy.profile(dom.bathy.profile>0);
    end    
    eta(j,:,:) = eta2temp;
    if model.phiForm==1
        phi(j,:,:) = phi2temp;
    else
        u(j,:,:)=u2temp;
        v(j,:,:)=v2temp;
    end
    if body.option==1
        for mm=1:Nbody
            RB(j,mm) = z_hat(j,2*nz+mm);
        end
    end
end
outkinematic=[];
if input.bdyassim.option==0
    % new space and fourier variables (needed if spaceinterpnum > 1)
    output.eta =eta;
    if strcmpi(options.OnlyEta,'No')
      if model.phiForm==1  
        output.phi =phi;
      else
        output.u   =u;
        output.v   =v;
      end
    else
        if model.phiForm==1
        output.phi =phi(end,:,:);
        else
        output.u   =u(end,:,:);
        output.v   =v(end,:,:);    
        end
    end
    output.time=time;
    output.X   =dom.X;
    output.Y   =dom.Y;
    output.break_nodes =dataBreak_nodes;
    output.break_crest =dataCrestBreak;
    output.break_speed =dataSpeedBreak;
    if body.option==1
        output.RB =RB;
        output.bodysavevar=bodysavevar;
    end
else
    if I<Npartition
        Ntn=length(time)-1;
    else
        Ntn=length(time);
    end
    Eta(ITERAssimPart:ITERAssimPart+Ntn-1,:,:)=eta(1:Ntn,:,:);
    if model.phiForm==1
    Phi(ITERAssimPart:ITERAssimPart+Ntn-1,:,:)=phi(1:Ntn,:,:);
    else
    U(ITERAssimPart:ITERAssimPart+Ntn-1,:,:)=u(1:Ntn,:,:); 
    V(ITERAssimPart:ITERAssimPart+Ntn-1,:,:)=v(1:Ntn,:,:);  
    end
    
    Time(ITERAssimPart:ITERAssimPart+Ntn-1)=time(1:Ntn);
    
    ITERAssimPart=ITERAssimPart+Ntn;
end

if checkInterior==1  
    global dteta dtphi timeIP phiIP etaIP
    if time(end)<par.IP.time(1)
        checkInterior=0;
    else
        if  Idstop==1   &&time(end)< par.IP.time(2)
            ind0=find(timeIP(2:end)==0,1,'first')+1;
            timeIP=timeIP(1:ind0-1,1);
            dteta=dteta(1:ind0-1,:,:);
            dtphi=dtphi(1:ind0-1,:,:);
            phiIP=phiIP(1:ind0-1,:,:);
            etaIP=etaIP(1:ind0-1,:,:);
        end
        
        try 
            if timeIP(end)<timeIP(end-1)
                ind0=find(timeIP(2:end)==0,1,'first')+1;
                dteta=dteta(1:ind0-1,:,:);
                dtphi=dtphi(1:ind0-1,:,:);
                etaIP=etaIP(1:ind0-1,:,:);
                phiIP=phiIP(1:ind0-1,:,:);
                timeIP=timeIP(1:ind0-1,1);
            end
        catch
        end
        outkinematic.dteta      =dteta;
        outkinematic.dtphi      =dtphi;
        outkinematic.eta        =etaIP;
        outkinematic.phi        =phiIP;
        outkinematic.time       =timeIP;
        if Npartition==1
            if options.mc.check==0
            save ('-v7.3',[Proj.workdir,Proj.savename,'_simul_Interior2D.mat'],'Proj','outkinematic','par','dom','Oprt','influx')
            else
            save ('-v7.3',[Proj.workdir,Proj.savename,'_simul_Interior2D_mc_',num2str(JJ),'.mat'],'Proj','outkinematic','par','dom','Oprt','influx')
            end
        else
           if options.mc.check==0
            save ('-v7.3',[Proj.workdir,Proj.savename,'_simul_Interior2D_part_',num2str(I),'.mat'],'Proj','outkinematic','par','dom','Oprt','influx')
           else
            save ('-v7.3',[Proj.workdir,Proj.savename,'_simul_Interior2D_mc_',num2str(JJ),'_part_',num2str(I),'.mat'],'Proj','outkinematic','par','dom','Oprt','influx')
           end
        
        end
        
    end
end


if Npartition==1
    output.comtime.relative=CompRel;
    output.comtime.ode=Timeodesolver;
    
    statusbarObj.setText(['saving data...']);
    
    if options.mc.check==0
        save ('-v7.3',[Proj.workdir,Proj.savename,'_simul','.mat'], ...
            'input','par','output','model','bath','dom','ivp','bdyassim','influx','Proj','timeSimul')
    else
        save ('-v7.3',[Proj.workdir,Proj.savename,'_simul_mc_',num2str(JJ),'.mat'], ...
            'input','par','output','model','bath','dom','ivp','bdyassim','influx','Proj','timeSimul')
    end
    statusbarObj.setText(['data saved']);
else
    if Idstop==1 || I==Npartition
        output.comtime.relative=CompRel;
        output.comtime.ode=Timeodesolver;
    end
    
    zinit_hat    = z_hat(end,:);
    
    if input.bdyassim.option==1
        if mod(I,floor(Npartition/Nsavedata))==0 || I==Npartition || Idstop==1
            output.eta =Eta;
            if model.phiForm==1
                if strcmpi(options.OnlyEta,'No')
                    output.phi =Phi;
                else
                    output.phi =phi(end,:,:);
                end
            else
                 if strcmpi(options.OnlyEta,'No')
                    output.u =U;output.v =V;
                else
                    output.u =u(end,:,:);
                    output.v =v(end,:,:);
                end
            end
            output.time=Time;
            output.X   =dom.X;
            output.Y   =dom.Y;
            output.break_nodes =dataBreak_nodes;
            output.break_crest =dataCrestBreak;
            output.break_speed =dataSpeedBreak;
            if body.option==1
                output.RB =RB;
                output.bodysavevar=bodysavevar;
            end
            
            statusbarObj.setText(['saving data...']);
            
            if Nsavedata==1
                if options.mc.check==0
                    save ('-v7.3',[Proj.workdir,Proj.savename,'_simul.mat'], ...
                        'input','par','output','model','bath','dom','ivp','bdyassim','influx','Proj','timeSimul')
                else
                    save ('-v7.3',[Proj.workdir,Proj.savename,'_simul_mc_',num2str(JJ),'.mat'], ...
                        'input','par','output','model','bath','dom','ivp','bdyassim','influx','Proj','timeSimul')
                end
                statusbarObj.setText(['data saved']);
            else
                if options.mc.check==0
                    save ('-v7.3',[Proj.workdir,Proj.savename,'_simul_part_',num2str(ITERAssimSave),'.mat'], ...
                        'input','par','output','model','bath','dom','ivp','bdyassim','influx','Proj','timeSimul')
                else
                    save ('-v7.3',[Proj.workdir,Proj.savename,'_simul_mc_',num2str(JJ),'_part_',num2str(ITERAssimSave),'.mat'], ...
                        'input','par','output','model','bath','dom','ivp','bdyassim','influx','Proj','timeSimul')
                    
                end
                statusbarObj.setText(['data part_',num2str(ITERAssimSave), ' saved']);
            end
            
            ITERAssimPart=1;ITERAssimSave=ITERAssimSave+1;
            
            if model.phiForm==1
            clearvars Eta Phi Time
            else
            clearvars Eta U V Time    
            end
            
            ITERbn=1; ITERdt=1;
            dataBreak_nodes=[];dataCrestBreak=[];
        end
        
    else
        statusbarObj.setText(['saving data...']);
        %  InputdataPostProcInternalFlows
        if options.mc.check==0
            save ('-v7.3',[Proj.workdir,Proj.savename,'_simul_part_',num2str(I),'.mat'], ...
                'input','par','output','model','bath','dom','ivp','bdyassim','influx','Proj','timeSimul')
            %     funC_savefast([Proj.workdir,Proj.savename,'_simul_part_',num2str(I),'.mat'], ...
            %         'input','par','output','model','bath','dom','ivp','bdyassim','influx','Proj','timeSimul')
            
        else
            save ('-v7.3',[Proj.workdir,Proj.savename,'_simul_mc_',num2str(JJ),'_part_',num2str(I),'.mat'], ...
                'input','par','output','model','bath','dom','ivp','bdyassim','influx','Proj','timeSimul')
            %     funC_savefast([Proj.workdir,Proj.savename,'_simul_mc_',num2str(JJ),'_part_',num2str(I),'.mat'], ...
            %         'input','par','output','model','bath','dom','ivp','bdyassim','influx','Proj','timeSimul')
            
        end
        statusbarObj.setText(['data part_',num2str(I), ' saved']);
    end
end
if strcmpi(options.OnlyEta,'Yes') && Npartition==1
    if options.mc.check==0
         if model.phiForm==1
          output.phi =phi;
         else
          output.u =u;
          output.v =v;
         end
    end
end

if input.bdyassim.option==1
    if strcmpi(options.OnlyEta,'Yes') && Npartition==1
        if options.mc.check==0
             if model.phiForm==1
             output.phi =Phi;
             else
              output.u =U;   
              output.v =V;   
             end
        end
    end
end
