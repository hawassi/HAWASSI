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
%%%%%%%%%% ODE-Solver                                             %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [output,outkinematic,log]=ODEsolver(h,preproc)

for JJ=1:length(preproc)
    %%Global parameters
    global ITERbn ITERdt iterProgBar flagbr dataBreak_nodes dataCrestBreak ...
        Idstop  iterdtcheck ITERbdt bodysavevar  itersaveBodyvar  dataSpeedBreak ...
        ITERbspdt
    iterdtcheck=1; iterProgBar=1; flagbr=0; ITERbn=1; ITERbspdt=1;
    dataBreak_nodes=[];dataCrestBreak=[]; ITERbdt=1; dataSpeedBreak=[];
    
    log=preproc(JJ).log;
    par=preproc(JJ).par;model=preproc(JJ).model;
    bath=preproc(JJ).bath;input=preproc(JJ).input;
    Proj=preproc(JJ).Proj;timeSimul=preproc(JJ).timeSimul;
    dom=preproc(JJ).dom;options=preproc(JJ).options;
    spatial=preproc(JJ).spatial;ivp=preproc(JJ).ivp;
    parBreak=preproc(JJ).parBreak;body=preproc(JJ).body;
    influx=preproc(JJ).influx;Oprt=preproc(JJ).Oprt;
    bdyassim=preproc(JJ).bdyassim;
    %     %%%adjust damping char func for user-defined fbl such that coef==1
    %     dom.fbl.char(dom.fbl.char==1/7)=dom.fbl.char(dom.fbl.char==1/7)/Oprt.Cpeak;
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%initialise fft method to find the most optimum method;
    fftw('planner','patient');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tintval=timeSimul.interval;
    checkInterior=options.interior.check;
    
    if checkInterior==1
        global iterInteriordt
        iterInteriordt=1;
        par.IP.time =options.interior.time;
        par.IP.check=1;
    else
        par.IP.check=0;
    end
    
    
    %%%initial condition
    if ivp.type==1 % non ivp
        
        if input.bdyassim.option==0 %influxing input
            zeta0_hat   = zeros(dom.Ny,dom.Nx);
            if model.phiForm==1
                zphi0_hat   = zeros(dom.Ny,dom.Nx);
                zinit_hat   = [reshape(zeta0_hat,1,[]),...
                    reshape(zphi0_hat,1,[])].';
            else
                zu0_hat   = zeros(dom.Ny,dom.Nx);
                zv0_hat   = zeros(dom.Ny,dom.Nx);
                zinit_hat   = [reshape(zeta0_hat,1,[]),...
                    reshape(zu0_hat,1,[]),reshape(zv0_hat,1,[])].';
            end
        else %%bdy assim input
            dtAssim=input.bdyassim.dt;
            tassimInit=input.bdyassim.tinit;tassimEnd=input.bdyassim.tend;
            
            bdyassimdata=input.bdyassim.assimdata;
            input.bdyassim.assimdata=[];%clear memory
            if tassimInit==tintval(1)
                Nx=bdyassimdata(2,3);
                Ny=bdyassimdata(3,3);
                Xdat=linspace(bdyassimdata(2,1),bdyassimdata(2,2),Nx);
                Ydat=linspace(bdyassimdata(3,1),bdyassimdata(3,2),Ny);
                bdyeta=bdyassimdata(4:3+Ny,:);
                EtaBdy=funBAssim_interpEtatoSimDom(Xdat,Ydat,bdyeta,dom);
                meandepthAssim=mean(mean(-dom.bathy.profile));
                
                if model.phiForm==1
                    if input.bdyassim.cb_phi==1
                        bdyassimdata_phi=input.bdyassim.assimdata_phi;
                        bdyphi=bdyassimdata_phi(4:3+Ny,:);
                        PhiBdy=funBAssim_interpEtatoSimDom(Xdat,Ydat,bdyphi,dom);
                        input.bdyassim.assimdata_phi=[];%clear memory
                    else
                        propdir=bdyassim.propdir;
                        PhiBdy=funOprt_PhifromEta(par.g,dom,model,meandepthAssim,EtaBdy,propdir);
                    end
                else
                    if input.bdyassim.cb_vel==1
                        bdyassimdata_u=input.bdyassim.assimdata_u;
                        bdyu=bdyassimdata_u(4:3+Ny,:);
                        uBdy=funBAssim_interpEtatoSimDom(Xdat,Ydat,bdyu,dom);
                        input.bdyassim.assimdata_u=[];%clear memory
                        
                        bdyassimdata_v=input.bdyassim.assimdata_v;
                        bdyv=bdyassimdata_v(4:3+Ny,:);
                        vBdy=funBAssim_interpEtatoSimDom(Xdat,Ydat,bdyv,dom);
                        input.bdyassim.assimdata_v=[];%clear memory
                    else
                        propdir=bdyassim.propdir;
                        PhiBdy=funOprt_PhifromEta(par.g,dom,model,meandepthAssim,EtaBdy,propdir);
                        gradphi=funOprt_grad2d(dom.Kx,dom.Ky,fft2(PhiBdy));
                        uBdy =gradphi.x;
                        vBdy =gradphi.y;
                    end
                    
                end
                
                
                EtaBdy=EtaBdy.*bdyassim.charupdate.*dom.cfSA;
                zeta0_hat   = fft2(EtaBdy);
                
                if model.phiForm==1
                    PhiBdy=PhiBdy.*bdyassim.charupdate.*dom.cfSA;
                    zphi0_hat   = fft2(PhiBdy);
                    zinit_hat   = [reshape(zeta0_hat,1,[]),...
                        reshape(zphi0_hat,1,[])].';
                else
                    uBdy=uBdy.*bdyassim.charupdate.*dom.cfSA;
                    vBdy=vBdy.*bdyassim.charupdate.*dom.cfSA;
                    
                    zu0_hat   = fft2(uBdy);
                    zv0_hat   = fft2(vBdy);
                    zinit_hat   = [reshape(zeta0_hat,1,[]),...
                        reshape(zu0_hat,1,[]),reshape(zv0_hat,1,[])].';
                end
                
            else
                zeta0_hat   = zeros(dom.Ny,dom.Nx);
                
                if model.phiForm==1
                    zphi0_hat   = zeros(dom.Ny,dom.Nx);
                    zinit_hat   = [reshape(zeta0_hat,1,[]),...
                        reshape(zphi0_hat,1,[])].';
                else
                    zu0_hat   = zeros(dom.Ny,dom.Nx);
                    zv0_hat   = zeros(dom.Ny,dom.Nx);
                    zinit_hat   = [reshape(zeta0_hat,1,[]),...
                        reshape(zu0_hat,1,[]),reshape(zv0_hat,1,[])].';
                end
                
            end
        end
    else      %Initial value problem;
        zeta0_hat   = fft2(ivp.eta);
        if model.phiForm==1
            zphi0_hat   = fft2(ivp.phi);
            zinit_hat   = [reshape(zeta0_hat,1,[]),...
                reshape(zphi0_hat,1,[])].';
        else
            zu0_hat   = fft2(ivp.u);
            zv0_hat   = fft2(ivp.v);
            zinit_hat   = [reshape(zeta0_hat,1,[]),...
                reshape(zu0_hat,1,[]),reshape(zv0_hat,1,[])].';
        end
    end
    
    
    
    
    bodysavevar=[];
    if body.option==1
        %%%preparation for saving variable in the ODE %%%%%%%%
        itersaveBodyvar=2;
        bodysavevar.betaz=zeros(influx.par.Nt,body.N);
        bodysavevar.betax=zeros(influx.par.Nt,body.N);
        bodysavevar.time=ones(influx.par.Nt,1).*influx.timesig(1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Nbody=body.N;
        zeta0_hat=fft2(funC_ifft2(zeta0_hat).*(1-body.par.chiB(end).char)+...
            body.par.shape.*body.par.chiB(end).char);
        xyz0                    = zeros(3*Nbody,1);
        xyz0(1:Nbody)           = body.xyz0(:,1);
        xyz0(Nbody+1:2*Nbody)   = body.xyz0(:,2);
        xyz0(2*Nbody+1:3*Nbody) = body.xyz0(:,3);
        orientationB            = zeros(3*Nbody,1);
        MomentumB               = zeros(6*Nbody,1);
        zinit_hat   = [reshape(zeta0_hat,1,[]),...
            reshape(zphi0_hat,1,[]),xyz0',orientationB',MomentumB'].';
    end
    
    %%preparing progress bar on the gui
    tic
    [jProgressBar,statusbarObj]=funGui_JavProgressBar(h);
    [jbStop]=funGui_Java_stopbutton(statusbarObj);
    set(jProgressBar,'Maximum',timeSimul.Nt, 'Value',0);
    jProgressBar.setStringPainted( true );
    ETA=funGui_remain_time(1,timeSimul.Nt);
    statusbarObj.setText(['estimating time remaining...']);
    if options.mc.check==1
        statusbarObj.setText(['Monte carlo simulation #',num2str(JJ),'>> estimating time remaining...']);
    end
    
    %%%Ode option
    diary on
    oo  = odeset('Events',@(time,z)funGui_eventStopfunction(time,z,jbStop),'RelTol',...
        par.ode.tol, 'AbsTol', par.ode.tol/10,'stats','on');% 1e-4 and 1e-6 standard
    
    
    ProgStatbar.jProgressBar=jProgressBar;
    ProgStatbar.statusbarObj=statusbarObj;
    
    if input.bdyassim.option==1
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
        
        Nsavedata=par.ode.Npartition;
        ITERAssimPart=1;ITERAssimSave=1;
        
    else
        Npartition=par.ode.Npartition;
    end
    
    
    
    if Npartition>1
        t_temp=tintval;
        Ntime_part=floor(length(t_temp)/Npartition);
    end
    
    Flag_timestep_failure=0;  %% flag time failure for Npartition>1
    
    DeltaTimesimul=0;
    disp('ODE statistic')
    cput=0;
    ITERbn=1; ITERdt=1;
    dataBreak_nodes=[];dataCrestBreak=[];
    
    
    for I=1:Npartition
        
        if Npartition>1
            flag       =0;
            if input.bdyassim.option==0
                ITERbn=1; ITERdt=1;
                dataBreak_nodes=[];dataCrestBreak=[];
            end
            
            if input.bdyassim.option==0
                if I<Npartition
                    tintval=t_temp((I-1)*Ntime_part+1:I*Ntime_part+1);
                end
                if I==Npartition
                    tintval=t_temp((I-1)*Ntime_part+1:end);
                end
            else
                tintval=tpartitionAss(I).time;
            end
            
            if I>1
                if input.bdyassim.option==1
                    disp(['Bdy. assimilation  #',num2str(I)])
                    
                    statusbarObj.setText(['Bdy. assimilation #',num2str(I)]);
                else
                    statusbarObj.setText(['Partition #',num2str(I)]);
                    disp(['Partition  #',num2str(I)])
                    
                end
            end
        end
        
        if  Idstop==1 && Npartition>1
            set(jProgressBar,'Maximum',timeSimul.Nt, 'Value',timeSimul.Nt);
            jProgressBar.setStringPainted( true );
            statusbarObj.setText('ODE calculation done.');
            Npartition=I-1;
            break;
        end
        
        
        
        if input.bdyassim.option==1 && I>1
            
            if I<=NAssimend
                Nx=bdyassimdata(2,3);
                Ny=bdyassimdata(3,3);
                
                if tinitSim<tassimInit
                    Iassim=I-1;
                else
                    Iassim=I;
                end
                
                Xdat=linspace(bdyassimdata(2,1),bdyassimdata(2,2),Nx);
                Ydat=linspace(bdyassimdata(3,1),bdyassimdata(3,2),Ny);
                bdyeta=bdyassimdata(4+(Iassim-1)*Ny:3+Iassim*Ny,:);
                EtaBdy=funBAssim_interpEtatoSimDom(Xdat,Ydat,bdyeta,dom);
                meandepthAssim=mean(mean(-dom.bathy.profile));
                
                if model.phiForm==1
                    if input.bdyassim.cb_phi==1
                        bdyphi=bdyassimdata_phi(4+(Iassim-1)*Ny:3+Iassim*Ny,:);
                        PhiBdy=funBAssim_interpEtatoSimDom(Xdat,Ydat,bdyphi,dom);
                    else
                        propdir=bdyassim.propdir;
                        PhiBdy=funOprt_PhifromEta(par.g,dom,model,meandepthAssim,EtaBdy,propdir);
                    end
                else
                    if input.bdyassim.cb_vel==1
                        bdyu=bdyassimdata_u(4+(Iassim-1)*Ny:3+Iassim*Ny,:);
                        bdyv=bdyassimdata_u(4+(Iassim-1)*Ny:3+Iassim*Ny,:);
                        uBdy=funBAssim_interpEtatoSimDom(Xdat,Ydat,bdyu,dom);
                        vBdy=funBAssim_interpEtatoSimDom(Xdat,Ydat,bdyv,dom);
                    else
                        propdir=bdyassim.propdir;
                        PhiBdy=funOprt_PhifromEta(par.g,dom,model,meandepthAssim,EtaBdy,propdir);
                        gradphi=funOprt_grad2d(dom.Kx,dom.Ky,fft2(PhiBdy));
                        uBdy =gradphi.x;
                        vBdy =gradphi.y;
                    end
                end
                
                EtaBdy=EtaBdy.*bdyassim.charupdate.*dom.cfSA;
                
                %               EtaBdy=EtaBdy.*bdyassim.chardata.*bdyassim.charupdate;
                %               PhiBdy=PhiBdy.*bdyassim.chardata.*bdyassim.charupdate;
                %
                z_hat=zinit_hat.';N=dom.Nx*dom.Ny;
                z1_hat=z_hat(1:N);     EtaLast_hat=reshape(z1_hat,dom.Ny,[]);
                EtaLast=EtaBdy+funC_ifft2(EtaLast_hat).*(1-bdyassim.charupdate);
                EtaLast_hat=fft2(EtaLast);
                
                if model.phiForm==1
                    PhiBdy=PhiBdy.*bdyassim.charupdate.*dom.cfSA;
                    z2_hat=z_hat(N+1:end); PhiLast_hat=reshape(z2_hat,dom.Ny,[]);
                    PhiLast=PhiBdy+funC_ifft2(PhiLast_hat).*(1-bdyassim.charupdate);
                    PhiLast_hat=fft2(PhiLast);
                else
                    uBdy=uBdy.*bdyassim.charupdate.*dom.cfSA;
                    vBdy=vBdy.*bdyassim.charupdate.*dom.cfSA;
                    
                    z2_hat=z_hat(N+1:2*N); uLast_hat=reshape(z2_hat,dom.Ny,[]);
                    z3_hat=z_hat(2*N+1:end); vLast_hat=reshape(z3_hat,dom.Ny,[]);
                    uLast=uBdy+funC_ifft2(uLast_hat).*(1-bdyassim.charupdate);
                    vLast=vBdy+funC_ifft2(vLast_hat).*(1-bdyassim.charupdate);
                    uLast_hat=fft2(uLast);vLast_hat=fft2(vLast);
                end
            else
                z_hat=zinit_hat.';N=dom.Nx*dom.Ny;
                z1_hat=z_hat(1:N);     EtaLast_hat=reshape(z1_hat,dom.Ny,[]);
                if model.phiForm==1
                    z2_hat=z_hat(N+1:end); PhiLast_hat=reshape(z2_hat,dom.Ny,[]);
                else
                    z2_hat=z_hat(N+1:2*N); uLast_hat=reshape(z2_hat,dom.Ny,[]);
                    z3_hat=z_hat(2*N+1:end);vLast_hat=reshape(z3_hat,dom.Ny,[]);
                end
            end
            zeta0_hat   = EtaLast_hat;
            if model.phiForm==1
                zphi0_hat   = PhiLast_hat;
                zinit_hat   = [reshape(zeta0_hat,1,[]),...
                    reshape(zphi0_hat,1,[])].';
            else
                zu0_hat   = uLast_hat;
                zv0_hat   = vLast_hat;
                zinit_hat   = [reshape(zeta0_hat,1,[]),...
                    reshape(zu0_hat,1,[]),reshape(zv0_hat,1,[])].';
            end
        end
        
        InternalPropPreparation2D;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% Evolution above FLAT/varying bottom (w/without wall)%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic;
        
        if body.option==0
            if  strcmpi(bath.name, 'Flat')
                if model.phiForm==1
                    if strcmpi(dom.wall.option,'No')
                        [time,z_hat] = ode45(@(time,z_hat)RHSFlat2d(time,z_hat,input,model,par,parBreak,dom,influx,ivp,Oprt,...
                            timeSimul,ProgStatbar),tintval,zinit_hat,oo);
                    else
                        [time,z_hat] = ode45(@(time,z_hat)RHSFlat2dWall(time,z_hat,input,model,par,parBreak,dom,influx,ivp,Oprt,...
                            timeSimul,ProgStatbar),tintval,zinit_hat,oo);
                    end
                else
                    if strcmpi(dom.wall.option,'No')
                        [time,z_hat] = ode45(@(time,z_hat)RHSFlat2duv(time,z_hat,input,model,par,parBreak,dom,influx,ivp,Oprt,...
                            timeSimul,ProgStatbar),tintval,zinit_hat,oo);
                    else
                        [time,z_hat] = ode45(@(time,z_hat)RHSFlat2dWalluv(time,z_hat,input,model,par,parBreak,dom,influx,ivp,Oprt,...
                            timeSimul,ProgStatbar),tintval,zinit_hat,oo);
                    end
                end
                
            elseif strcmpi(bath.name,'Shore')
                %%update initial condition
                %                     [zeta_b0]   = InitValEtaRunUp2D(dom.bathy.profile,dom.bathy.HminShore);
                %                     zeta0       = funC_ifft2(zeta0_hat);
                %                     zeta0_hat   = fft2(zeta0+zeta_b0);
                
                %                     if ~strcmp(influx.par.category,'Shallow water')
                %                         [zphi_b0]   = InitValPhiRunUp2D(dom.bathy.profile,dom.bathy.HminShore,dom.dx,dom.dy);
                %                         zphi0       = funC_ifft2(zphi0_hat);
                %                         zphi0_hat   = fft2(zphi0+zphi_b0);
                %                     end
                
                
                if model.phiForm==1
                    [time,z_hat] = ode45(@(time,z_hat)RHSShore2d_init0(time,z_hat,input,model,par,bath,parBreak,dom,influx,ivp,Oprt,...
                        timeSimul,ProgStatbar),tintval,zinit_hat,oo);
                else
                    
                    [time,z_hat] = ode45(@(time,z_hat)RHSShore2duv_init0(time,z_hat,input,model,par,bath,parBreak,dom,influx,ivp,Oprt,...
                        timeSimul,ProgStatbar),tintval,zinit_hat,oo);
                end
                %                     [time,z_hat] = ode45(@(time,z_hat)RHSShore2d(time,z_hat,input,model,bath,par,parBreak,dom,influx,ivp,Oprt,...
                %                                  timeSimul,ProgStatbar),tintval,zinit_hat,oo);
                
                
            else %%varying bottom
                if strcmpi(dom.wall.option,'No')
                    if model.phiForm==1
                        [time,z_hat] = ode45(@(time,z_hat)RHSBathy2d(time,z_hat,input,model,par,parBreak,dom,influx,ivp,Oprt,...
                            timeSimul,ProgStatbar),tintval,zinit_hat,oo);
                    else
                        [time,z_hat] = ode45(@(time,z_hat)RHSBathy2duv(time,z_hat,input,model,par,parBreak,dom,influx,ivp,Oprt,...
                            timeSimul,ProgStatbar),tintval,zinit_hat,oo);
                    end
                else
                    if model.phiForm==1
                        [time,z_hat] = ode45(@(time,z_hat)RHSBathy2dWall(time,z_hat,input,model,par,parBreak,dom,influx,ivp,Oprt,...
                            timeSimul,ProgStatbar),tintval,zinit_hat,oo);
                    else
                        [time,z_hat] = ode45(@(time,z_hat)RHSBathy2dWalluv(time,z_hat,input,model,par,parBreak,dom,influx,ivp,Oprt,...
                            timeSimul,ProgStatbar),tintval,zinit_hat,oo);
                    end
                end
            end
        else
            [time,z_hat] = ode45(@(time,z_hat)RHSBody2d(time,z_hat,input,model,par,dom,influx,ivp,body,Oprt,...
                timeSimul,ProgStatbar),tintval,zinit_hat,oo);
        end
        cputf=toc;
        cput=cput+cputf;
        DeltaTimesimul=DeltaTimesimul+(time(end)-time(1));
        
        ODEoutputprocessor;
    end
    jbStop.setVisible(0);
    jProgressBar.setVisible(0);
    
    if input.bdyassim.option==0
        if Npartition>1
            if options.partition.combine==1
                Combine_output_part2d;
            end
        end
    else
        if ITERAssimSave-1>1
            if options.partition.combine==1
                Combine_output_part2d;
            end
        end
    end
    
    N0=length(log.string);N=N0;
    log.string{N+1}=['Timeodesolver= ',num2str(Timeodesolver)];
    log.string{N+2}=['CompRel= ',num2str(CompRel)];
    N=N+2;
    if Flag_timestep_failure==1
        log.string{N+1}='Time stepping failure!';
        N=N+1;
    end
    log.string{N+1}='==========================================================';
    disp('--------------------------------------------------------------------');
    disp(log.string{N0+1});
    disp(log.string{N0+2});
    disp(log.string{N0+3});
    if Flag_timestep_failure==1
        disp(log.string{N0+4});
    end
    diary off;
end
end
