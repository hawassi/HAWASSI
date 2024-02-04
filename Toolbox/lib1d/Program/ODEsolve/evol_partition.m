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
%%%%%%%%%%%%%%%%%  Evolution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Global Parameter
global ITERdt flag dataCrestBreak dataBreak_nodes flagWarnIVP ITERbn iterProgBar ...
    Idstop ITERbdt tprev shipsavevar itersavevar Sxprev  EvL EvR
flagWarnIVP=1; ITERbn=1; ITERbdt=1; ITERdt=1; iterProgBar=1;
flag       =0; dataBreak_nodes=[];dataCrestBreak=[]; tprev=0;

               


cutfrac = par.cutfrac; k= par.k; x = par.x; k_p  = influx.k_p; dtout = par.dt;
depth_inf  = influx.depth;
g          = par.g;
par.Nt     = length(influx.gen.timesig);
Npartition = par.ode_Npartition;
checkInterior=options.interior.check;


if checkInterior==1
    global iterInteriordt
    iterInteriordt=1; %%variable flag for partition ODE
    par.IP.time =options.interior.time;
    par.IP.check=1;
else
    par.IP.check=0;
end


if IVP.type==1 % non ivp
    if input.assim.check==0 % influxing
        zeta0_hat   = zeros(length(k),1);
        zu0_hat     = zeros(length(k),1);
        zini        = [zeta0_hat,zu0_hat];
    else %% assimilation
        tassimInit=input.assim.data.eta(2,1);tassimEnd=input.assim.data.eta(end,1);
        bdyassimdata=input.assim.data.eta;
        input.assim.data.eta=[];%clear memory
        if tassimInit==tintval(1)
            xdat=bdyassimdata(1,2:end);
            bdyeta=bdyassimdata(2,2:end);
            EtaBdy=interp1(xdat,bdyeta,par.x,'spline');
            
            meandepthAssim=mean(-par.bathy);
            if input.assim.data.checkbox_velocity==1
                bdyassimdata_phi=input.assim.data.velocity;
                bdyu=bdyassimdata_phi(2,2:end);
                UBdy=interp1(xdat,bdyu,par.x,'spline');
                input.bdyassim.assimdata_phi=[];%clear memory
            else
                propdir=bdyassim.propdir;
                phihat=fft(EtaBdy)./(1i.*OmExact(par.k.',meandepthAssim,[]).' );
                phihat(abs(phihat)==Inf)=0;
                UBdy=g*Ifft(1i.*par.k.'.*phihat);
                if propdir==2
                    UBdy=-UBdy;
                end
            end
            EtaBdy=EtaBdy.*bdyassim.charupdate.*par.cfSA;
            UBdy=UBdy.*bdyassim.charupdate.*par.cfSA;
            
            zeta0_hat   = fft(EtaBdy);
            zu0_hat   = fft(UBdy);
            zini        = [zeta0_hat,zu0_hat];
        else
            zeta0_hat   = zeros(length(k),1);
            zu0_hat     = zeros(length(k),1);
            zini        = [zeta0_hat,zu0_hat];
        end
    end
     
else      %Initial value problem;
    zeta0_hat   = IVP.zeta0_hat;
    zu0_hat     = IVP.zu0_hat;
    zini        = [zeta0_hat,zu0_hat];
end

if shippar.check== 1
    %%%% global variables
   
    EvL=0; EvR=0;
    shipsavevar=[];
    itersavevar=2;
    Nt=floor((par.t_end-par.t_init)./par.dt)+1;
    shipsavevar.zVel=zeros(Nt,shippar.Nship);
    shipsavevar.xVel=zeros(Nt,shippar.Nship);
    shipsavevar.thetaVel=zeros(Nt,shippar.Nship);
    shipsavevar.betax=zeros(Nt,shippar.Nship);
    shipsavevar.betaz=zeros(Nt,shippar.Nship);
    shipsavevar.betatheta=zeros(Nt,shippar.Nship);
    shipsavevar.dxi_Kdiffx=zeros(Nt,shippar.Nship);
    shipsavevar.dxi_Kdiffz=zeros(Nt,shippar.Nship);
    shipsavevar.dxi_Kdifftheta=zeros(Nt,shippar.Nship);
    shipsavevar.dxi_Kradx=zeros(Nt,shippar.Nship);
    shipsavevar.dxi_Kradz=zeros(Nt,shippar.Nship);
    shipsavevar.dxi_Kradtheta=zeros(Nt,shippar.Nship);
    
    shipsavevar.phiwave=zeros(Nt,length(x));
    shipsavevar.phirad=zeros(Nt,length(x));
    shipsavevar.addMass=zeros(Nt,9);
    shipsavevar.dampCoef=zeros(Nt,9);
    shipsavevar.time=ones(Nt,1).*par.t_init;
    
    Nship       = shippar.Nship;
     if strcmp(bath.type,'BR') || strcmp(bath.type,'UR')
     etaI0       =shippar.init.zeta0;
     etaI0((etaI0-par.bathy)<=par.interp.H_minShore)=0;
      else
     etaI0       =shippar.init.zeta0;
     end
   
    zeta0_hat   = fft(etaI0);
    zphi0_hat   = fft(shippar.init.phi0);%fft(Ifft(phihat).*ship.form.chi(:,end));
%     figure;
%     subplot(2,1,1)
%     plot(etaI0)
%     subplot(2,1,2)
%     plot(shippar.init.phi0)
%     
    shipsavevar.phiwave(1,:)=shippar.init.phi0;
    
    Sxprev=shippar.form.xShip;
    RB0         = zeros(6*Nship,1);
  
    try
    RB0(1:Nship)=0;%x(closest(x,str2num(cell2mat(shippar.data(:,5))))).'-shippar.form.xShip(:,2);
    RB0(2*Nship+1:3*Nship)=str2num(cell2mat(shippar.data(:,6)));   
    catch
    RB0(1:Nship)=0;%x(closest(x,cell2mat(shippar.data(:,5)))).'-shippar.form.xShip(:,2);
    RB0(2*Nship+1:3*Nship)=  cell2mat(shippar.data(:,6));
    end
    if IVP.type==4
    RB0(1:6*Nship,1)=IVP.data(1:6*Nship,4);   
    end
   
    Extra0      = zeros(6*Nship,1);
    zini        = [zeta0_hat;zphi0_hat;RB0;Extra0];
end

tic
JavProgressBar;
[jbStop]=Java_stopbutton(statusbarObj);
set(jProgressBar,'Maximum',par.Nt, 'Value',0);
jProgressBar.setStringPainted( true );
ETA=remain_time(1,par.Nt);
statusbarObj.setText(['estimating time remaining...']);
oo  = odeset('Events',@(time,z)eventStopfunction(time,z,jbStop),'RelTol',...
    par.odetol, 'AbsTol', par.odetol/10);% 1e-4 and 1e-6 standard


aal  = fun_alias(k,cutfrac); %wavenumbercut (anti-alias)

par.Nt=length(tintval);


ProgStatbar.jProgressBar=jProgressBar;
ProgStatbar.statusbarObj=statusbarObj;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%initialise fft method to find the most optimum method;
fftw('planner','patient');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
if input.assim.check==1
    dtAssim=input.assim.data.tstep;
    Npartition=floor((tassimEnd-tassimInit)/dtAssim);
    tinitSim=tintval(1);
    if tintval(end)>tassimEnd
        Npartition=Npartition+1;
    end
    
    if tintval(1)<tassimInit
        Npartition=Npartition+1;
    end
    
    if tintval(end)>tassimEnd
        NAssimend=Npartition-1;
    else
        NAssimend=Npartition;
    end
    
    tpartitionAss=funBAssim_tpartition(tintval,tassimInit,tassimEnd,dtAssim,Npartition);
    Nsavedata=par.ode_Npartition;
    ITERAssimPart=1;ITERAssimSave=1;  
else
    Npartition=par.ode_Npartition;
end

if Npartition>1
    t_temp=tintval;
    Ntime_part=floor(length(t_temp)/Npartition);
end
Flag_timestep_failure=0;  %% flag time failure for Npartition>1
tic;
% figure;


for I=1:Npartition
    if Npartition>1
        ITERbn=1; ITERdt=1; iterProgBar=1;
        flag       =0; dataBreak_nodes=[];dataCrestBreak=[];
        if input.assim.check==0
            if I<Npartition
                tintval=t_temp((I-1)*Ntime_part+1:I*Ntime_part+1);
            else
                tintval=t_temp((I-1)*Ntime_part+1:end);
            end
        else
            tintval=tpartitionAss(I).time;
        end
    end
    
    if I>1
        if input.assim.check==1
            %   disp(['Bdy. assimilation  #',num2str(I)])
            
            statusbarObj.setText(['Bdy. assimilation #',num2str(I)]);
        else
            statusbarObj.setText(['Partition #',num2str(I)]);
            %    disp(['Partition  #',num2str(I)])
            
        end
    else
        if Npartition>1
            if input.assim.check==1
                %    disp(['Bdy. assimilation  #',num2str(I)])
                
                statusbarObj.setText(['Bdy. assimilation #',num2str(I)]);
            else
                statusbarObj.setText(['Partition #',num2str(I)]);
                %     disp(['Partition  #',num2str(I)])
                
            end
        end
    end
    
    
    if  Idstop==1 && Npartition>1
        set(jProgressBar,'Maximum',par.Nt, 'Value',par.Nt);
        jProgressBar.setStringPainted( true );
        statusbarObj.setText('ODE calculation done.');
        Timeodesolver =toc;
        CompRel = 100*Timeodesolver/(par.t_end - par.t_init);
        
        Npartition=I-1;
        break;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if input.assim.check==1 && I>1
         if I<=NAssimend
            if tinitSim<tassimInit
                Iassim=I-1;
            else
                Iassim=I;
            end
            
            bdyeta=bdyassimdata(1+Iassim,2:end);
            EtaBdy=interp1(xdat,bdyeta,par.x,'spline');
            meandepthAssim=mean(mean(-par.bathy));
            
            if input.assim.data.checkbox_velocity==1
                bdyu=bdyassimdata_phi(1+Iassim,2:end);
                UBdy=interp1(xdat,bdyu,par.x,'spline');
            else
                phihat=fft(EtaBdy)./(1i.*OmExact(par.k.',meandepthAssim,[]).' );
                phihat(abs(phihat)==Inf)=0;
                UBdy=g*Ifft(1i.*par.k.'.*phihat);
                if propdir==2
                    UBdy=-UBdy;
                end
            end
            
            EtaBdy=EtaBdy.*bdyassim.charupdate.*par.cfSA;
            UBdy=UBdy.*bdyassim.charupdate.*par.cfSA;
            %               EtaBdy=EtaBdy.*bdyassim.chardata.*bdyassim.charupdate;
            %               PhiBdy=PhiBdy.*bdyassim.chardata.*bdyassim.charupdate;
            %
            z_hat=zini;N=par.Nx;
            EtaLast_hat=z_hat(1:N);    
            ULast_hat=z_hat(N+1:end);
            EtaLast=EtaBdy+Ifft(EtaLast_hat).*(1-bdyassim.charupdate);
            ULast=UBdy+Ifft(ULast_hat).*(1-bdyassim.charupdate);
            EtaLast_hat=fft(EtaLast); ULast_hat=fft(ULast);
%             plot(x,EtaLast,'r',x,ULast,'b');
%             tnow=tintval;
%             pause(0.001);
        else
            z_hat=zini;N=par.Nx;
            EtaLast_hat=z_hat(1:N);
            ULast_hat=z_hat(N+1:end);
        end
        zeta0_hat   = EtaLast_hat;
        zu0_hat     = ULast_hat;
        zini        = [zeta0_hat,zu0_hat];
   end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    InternalPropPreparation;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Evolution above FLAT bottom (w/without wall) %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if shippar.check== 0
        if strcmp(bath.type,'F')
           if strcmp(bath.wall.check,'No')


                [time,z] = ode45(@(time,z)RHSFlat(time,z,model,par,influx,IVP,Oprt,...
                    ProgStatbar),tintval,zini,oo);
           else
              if par.wall.method==1
                 [time,z] = ode45(@(time,z)RHSFlat(time,z,model,par,influx,IVP,Oprt,...
                    ProgStatbar),tintval,zini,oo); 
              else
                if bath.wall.type==1
                    [time,z] = ode45(@(time,z)RHSFlatWall(time,z,model,par,influx,IVP,Oprt,...
                        ProgStatbar),tintval,zini,oo);
                else
                    [time,z] = ode45(@(time,z)RHSFlatWallFreq(time,z,model,par,influx,IVP,Oprt,...
                        ProgStatbar),tintval,zini,oo);
                end
              end
                
           end
            
        end
        
        
        if strcmp(bath.type,'B') || strcmp(bath.type,'U')
            if strcmp(bath.wall.check,'No')
                [time,z] = ode45(@(time,z)RHSBathy(time,z,model,par,influx,IVP,Oprt,...
                    ProgStatbar),tintval,zini,oo);
            else
                if par.wall.method==1
                  [time,z] = ode45(@(time,z)RHSBathy(time,z,model,par,influx,IVP,Oprt,...
                    ProgStatbar),tintval,zini,oo);   
                else
                    if bath.wall.type==1
                        [time,z] = ode45(@(time,z)RHSBathyWall(time,z,model,par,influx,IVP,Oprt,...
                            ProgStatbar),tintval,zini,oo);
                    else
                        [time,z] = ode45(@(time,z)RHSBathyWallFreq(time,z,model,par,influx,IVP,Oprt,...
                            ProgStatbar),tintval,zini,oo);
                    end
                end
            end
        end
        
        if strcmp(bath.type,'BR') || strcmp(bath.type,'UR')
            if I==1
             
%                 [zeta_b0]   = InitValEtaRunUp(par.bathy,par.interp.H_minShore);
%                 zeta0       =  Ifft(zeta0_hat);
%                 zeta0_hat       = fft(zeta0+zeta_b0);
               

                
                zini        = [zeta0_hat;zu0_hat];
                
            end
        
%             [time,z] = ode45(@(time,z)RHSShore(time,z,model,par,influx,IVP,Oprt,...
%                 ProgStatbar),tintval,zini,oo);
            
            [time,z] = ode45(@(time,z)RHSShore_etainitnol1_smoothchar(time,z,model,par,influx,IVP,Oprt,...
                ProgStatbar),tintval,zini,oo);
            
%              [time,z] = ode45(@(time,z)RHSShore_etainitnol1_phi(time,z,model,par,influx,IVP,Oprt,...
%                 ProgStatbar),tintval,zini,oo);

%             [time,z] = ode45(@(time,z)RHSShore180801(time,z,model,par,influx,IVP,Oprt,...
%                 ProgStatbar),tintval,zini,oo);

        end
        
    else
        if strcmp(bath.type,'F')
            if strcmp(bath.wall.check,'No')
                [time,z] = ode45(@(time,z)RHSFlatShip2etavirt(time,z,model,par,influx,IVP,Oprt,...
                    shippar,ProgStatbar),tintval,zini,oo);
             else
                if par.wall.method==1
                  [time,z] = ode45(@(time,z)RHSFlatShip2etavirt(time,z,model,par,influx,IVP,Oprt,...
                    shippar,ProgStatbar),tintval,zini,oo);  
                else
                   [time,z] = ode45(@(time,z)RHSFlatShipWall1(time,z,model,par,influx,IVP,Oprt,...
                    shippar,ProgStatbar),tintval,zini,oo); 
                end
            end
                
        end
        
        if strcmp(bath.type,'B') || strcmp(bath.type,'U')
            if strcmp(bath.wall.check,'No')
                [time,z] = ode45(@(time,z)RHSBathyShip2etavirt(time,z,model,par,influx,IVP,Oprt,...
                    shippar,ProgStatbar),tintval,zini,oo);
            else
                
            end
        end
        
         if strcmp(bath.type,'BR') || strcmp(bath.type,'UR')
            if I==1
%                 [zeta_b0]   = InitValEtaRunUp(par.bathy,par.interp.H_minShore);
%                 zeta0       =  Ifft(zeta0_hat)*0;
%                 zeta0_hat       = fft(zeta0+zeta_b0);
                zini        = [zeta0_hat;zphi0_hat;RB0;Extra0];
            end
            
            
            %if strcmp(bath.wall.check,'No')
            [time,z] = ode45(@(time,z)RHSShoreShip_etavirt(time,z,model,par,influx,IVP,Oprt,...
                shippar,ProgStatbar),tintval,zini,oo);
%             %else 
            %end
        end
        
    end
    

  
    
    if Npartition>1
        if Idstop==0
            if time(end)<tintval(end)  %time step failure
                Idstop=1;
                Flag_timestep_failure=1;
            end
        end
    end
    
    
    if Npartition==1
        set(jProgressBar,'Maximum',par.Nt, 'Value',par.Nt);
        jProgressBar.setStringPainted( true );
        statusbarObj.setText('ODE calculation done.');
        statusbarObj.setText(['saving data...']);
        Timeodesolver =toc;
        CompRel = 100*Timeodesolver/(par.t_end - par.t_init);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Producing output of ode-calculation%
    Nend        = length(time);
    if shippar.check== 1
        nz          = (length(z(1,:))-2*6*Nship)/2;
        RB          = zeros(Nend,6*Nship);
        Ext         = zeros(Nend,6*Nship);
    else
        nz          = length(z(1,:))/2;
    end
    
    eta       = zeros(Nend,nz);
    u         = zeros(Nend,nz);
    for j = 1 : length(time)
        eta2hat  = z(j,1:nz);
        eta(j,:) = Ifft(eta2hat);
       
        if strcmp(bath.type,'BR') || strcmp(bath.type,'UR')
         eta(j,:) =eta(j,:)+par.bathyplus.';   
        end
        if strcmp(par.wall.presence,'Yes') && par.wall.method==1
            if par.wall.type==1
                eta(j,closest(par.x,par.wall.position))=(1+par.wall.refl_Coef)*eta(j,closest(par.x,par.wall.position));
            else
                Rho_k_Eta=Ifft(fft(eta(j,:)).*(1+par.wall.refl_Coef_k'));
                eta(j,closest(par.x,par.wall.position))=Rho_k_Eta(closest(par.x,par.wall.position));
            end 
        end
        if  strcmp(options.OnlyEta,'No')
            if shippar.check== 0
                u2hat  = z(j,nz+1:2*nz);
            else
                phi_hat_j=z(j,nz+1:2*nz);
                u2hat=1i.*par.k'.*phi_hat_j;
            end
            u(j,:)= Ifft(u2hat);
            if strcmp(bath.wall.check,'Yes')
                if par.wall.type==1
                    u(j,closest(par.x,par.wall.position))=(1+par.wall.refl_Coef)*u(j,closest(par.x,par.wall.position));
                else
                    Rho_k_U=Ifft(fft(u(j,:)).*(1+par.wall.refl_Coef_k'));
                    u(j,closest(par.x,par.wall.position))=Rho_k_U(closest(par.x,par.wall.position));
                end
            end
        end
        
        if shippar.check== 1
            for mm=1:6*Nship
                RB(j,mm) = z(j,2*nz+mm);
                Ext(j,mm)= z(j,2*nz+6*Nship+mm);
            end
        end
    end

  
    dataCrestBreak;dataBreak_nodes;
    
    dx      = par.dx;
    x       = par.x;
    
    k     =   freqspace(x);
    if input.assim.check==0
        output.eta =eta;
        output.u   =u;
        output.time=time;
        output.x   =x;
        output.break_nodes =dataBreak_nodes;
        output.break_crest =dataCrestBreak;
        if shippar.check== 1
            output.RB=RB;
            output.Ext=Ext;
            output.shipsavevar=shipsavevar;
            assignin('base','shipsavevar',shipsavevar)
        end
    else
        if I<Npartition
            Ntn=length(time)-1;
        else
            Ntn=length(time);
        end
        Eta(ITERAssimPart:ITERAssimPart+Ntn-1,:)=eta(1:Ntn,:);
        U(ITERAssimPart:ITERAssimPart+Ntn-1,:)=u(1:Ntn,:);
        Time(ITERAssimPart:ITERAssimPart+Ntn-1)=time(1:Ntn);
        
        ITERAssimPart=ITERAssimPart+Ntn;
    end

  
    if Npartition==1
        statusbarObj.setText(['saving data...']);
        
        InputdataPostProcInternalFlows;
        
        save ('-v7.3',[Proj.Dir,Proj.savename,'_simul','.mat'], ...
            'output','model','bath','par','influx','Proj','shippar','Oprt')
        
        statusbarObj.setText(['data saved']);
    else
        
        
        zini        = z(end,:);
        
        if input.assim.check==1

            if mod(I,floor(Npartition/Nsavedata))==0 || I==Npartition || Idstop==1
                output.eta =Eta;
                if strcmpi(options.OnlyEta,'No')
                    output.u =U;
                end
                output.time=Time;
                output.x   =x;
                output.break_nodes =dataBreak_nodes;
                output.break_crest =dataCrestBreak;
                if shippar.check== 1
                    output.RB=RB;
                    output.Ext=Ext;
                    output.shipsavevar=shipsavevar;
                end
                statusbarObj.setText(['saving data...']);
                InputdataPostProcInternalFlows;
                
                if Nsavedata==1
                    save ('-v7.3',[Proj.Dir,Proj.savename,'_simul.mat'], ...
                        'output','model','bath','par','influx','Proj','shippar','Oprt')
                    statusbarObj.setText(['data saved']);
                else
                    save ('-v7.3',[Proj.Dir,Proj.savename,'_simul_part_',num2str(I),'.mat'], ...
                        'output','model','bath','par','bdyassim','influx','Proj','shippar','Oprt')
                    statusbarObj.setText(['data part_',num2str(I), ' saved']);
                end
                
                ITERAssimPart=1;ITERAssimSave=ITERAssimSave+1;
                
                clearvars Eta Phi Time
                
                ITERbn=1; ITERdt=1;
                dataBreak_nodes=[];dataCrestBreak=[];         
            end
            
        else
            statusbarObj.setText(['saving data...']);
            InputdataPostProcInternalFlows;
            save ('-v7.3',[Proj.Dir,Proj.savename,'_simul_part_',num2str(I),'.mat'], ...
                'output','model','bath','par','influx','Proj','shippar','Oprt')
            statusbarObj.setText(['data part_',num2str(I), ' saved']);
        end
    end
   
  
end

if Npartition>1
    set(jProgressBar,'Maximum',par.Nt, 'Value',par.Nt);
    jProgressBar.setStringPainted( true );
    statusbarObj.setText('ODE calculation done.');
    statusbarObj.setText(['saving data...']);
    Timeodesolver =toc;
    CompRel = 100*Timeodesolver/(par.t_end - par.t_init);
end
jbStop.setVisible(0);
jProgressBar.setVisible(0);


if input.assim.check==0
    if Npartition>1
        if options.partition.combine==1
            combines_part_output;
        end
    end
else
    if ITERAssimSave-1>1
        if options.partition.combine==1
            combines_part_output;
        end
    end
end

  