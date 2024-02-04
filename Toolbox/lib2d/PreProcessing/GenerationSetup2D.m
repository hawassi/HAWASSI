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
%%%%%%%%%% Generation Setup                                      %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OmExct          = str2func('funOprt_OmExact');
UgExct          = str2func('funOprt_UgExact');
omAdd           = model.OmFun;
ugAdd           = model.UgFun;
g               = par.g;

if strcmpi(input.wave.option,'Yes')
[input.wave,model.influx,spatial.influx]=funG_adjust_influx_param_from_gui(input.wave,dom);
assignin('base','spatinfl',spatial.influx)
lambda_p        =zeros(input.wave.N,1);
k_p             =zeros(input.wave.N,1);
Tp              =zeros(input.wave.N,1);
Tm01            =zeros(input.wave.N,1);
Hs              =zeros(input.wave.N,1);
nupeak          =zeros(input.wave.N,1);
meandepthinflux =zeros(input.wave.N,1);
wavecateg       =cell(input.wave.N);
ChiAdj          =1;
Cg_p            =zeros(input.wave.N,1);


for I=1:input.wave.N
    inflXY=spatial.influx.line(I).xy;
    depthinfluxline=-dom.bathy.profile(inflXY(:,3));
    meandepthinflux(I)     = mean(depthinfluxline);
   
    if input.wave.type(I).flag==2 || input.wave.type(I).flag ==3
        Tp(I)       = input.wave.Tp(I);            % [s] Peak period for JS or period for harmonic
        Hs(I)       = input.wave.Hs(I);            % = significant wave height [m] for JS
        ampl        = Hs(I)/2;
        dtsig       = input.wave.dt(I);%Tp/30;        % choice of accuracy input signal
        t_init      = input.wave.t_init(I);
        t_end       = input.wave.t_end(I);
        timesig     = [t_init:dtsig:t_end];
        Ntsig       = length(timesig);
        if  mod(Ntsig,2)~=0  %Ntsig must even for input JonsWap
            timesig     = [t_init:dtsig:t_end-dtsig];
            Ntsig       = length(timesig);
        end
        omsig       = funC_freqspace(timesig);
        halfomsig   = omsig(1:floor(end/2))';
    end
    
    %%%%%%%%%%%%%% Setup generation signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if  input.wave.type(I).flag  ==  2 %harmonic
        theta0=deg2rad(input.wave.direction(I));
        [insigI,influx.Spect(I).prop]=funG_Harmonic_2D_signal(g,ampl,Tp(I),dom,timesig,depthinfluxline,spatial.influx,theta0,I);
    elseif input.wave.type(I).flag  ==  3 %JONSWAP
        [insigI,influx.Spect(I).prop]=funG_JONSWAP_2D_signal(g,dom,input.wave,timesig,halfomsig,spatial.influx,depthinfluxline,I);
    elseif input.wave.type(I).flag  ==  4 %User defined (signal)
        [insigI,timesig,omsig,halfomsig]=funG_UserDefinedSignalInterp2(input.wave,dom,spatial.influx,I);
    elseif input.wave.type(I).flag  ==  5 %User defined (signal)
        [insigI,timesig,omsig,halfomsig,influx.Spect(I).prop]=funG_UserDefinedSpectrum(g,input.wave,dom,spatial.influx,depthinfluxline,I); 
    end   
 
     
                
    %%Saving influx data for input type 2,3,4&5
       
       InfluxSignal=zeros(length(timesig)+2,length(inflXY(:,1))+1);
       InfluxSignal(1,2:end)=inflXY(:,1);
       InfluxSignal(2,2:end)=inflXY(:,2);
       if strcmpi(spatial.influx.line(I).Orientation,'Horizontal')
       Ind1=funC_closest(dom.X,min(inflXY(:,1)));
       Ind2=funC_closest(dom.X,max(inflXY(:,1)));
       else
       Ind1=funC_closest(dom.Y,min(inflXY(:,2)));
       Ind2=funC_closest(dom.Y,max(inflXY(:,2)));
       end
       InfluxSignal(3:end,:)=[timesig' insigI(:,Ind1:Ind2)];
       
       [OmS,S_hat]=funSP_variance_density_spectrum1d(timesig,insigI(:,floor((Ind1+Ind2)/2)));
       
       VarDensSpect1D=zeros(length(OmS),2);
       VarDensSpect1D(:,1)=abs(OmS);VarDensSpect1D(:,2)=S_hat;
       
       if input.wave.type(I).flag  ==  3 || input.wave.type(I).flag  ==  5
        halfomsig=influx.Spect(I).prop.halfomsig;  
        mdir_rad=influx.Spect(I).prop.mDir;
        cos2pdf=influx.Spect(I).prop.cos2pdf;
        varDensSpect=influx.Spect(I).prop.varDensSpect;
        stdTheta=rad2deg(influx.Spect(I).prop.StdTheta);
        theta0  = linspace(-pi/2+mdir_rad,pi/2+mdir_rad, length(cos2pdf));
        Spect=varDensSpect'*cos2pdf;
        VarDensSpect2D=zeros(length(halfomsig)+1,length(theta0)+1);
        VarDensSpect2D(2:end,1)=halfomsig;
        VarDensSpect2D(1,2:end)=theta0;
        VarDensSpect2D(2:end,2:end)=Spect;
       end
       InflDir=[Proj.workdir,'\InfluxData\'];
     
       if ~isdir(InflDir), mkdir(InflDir); end
       if options.mc.check==0
       save ('-v7.3',[InflDir,'InfluxSignal_',num2str(I),'.mat'],'InfluxSignal')
       save ('-v7.3',[InflDir,'VarDensSpect1D_',num2str(I),'.mat'],'VarDensSpect1D')
        if input.wave.type(I).flag  ==  3 || input.wave.type(I).flag  ==  5
        save ('-v7.3',[InflDir,'VarDensSpect2D_',num2str(I),'.mat'],'VarDensSpect2D')
        end
       else
       save ('-v7.3',[InflDir,'InfluxSignal_',num2str(I),'_mc',num2str(JJ),'.mat'],'InfluxSignal')
       save ('-v7.3',[InflDir,'VarDensSpect1D_',num2str(I),'_mc',num2str(JJ),'.mat'],'VarDensSpect1D')
        if input.wave.type(I).flag  ==  3 || input.wave.type(I).flag  ==  5
        save ('-v7.3',[InflDir,'VarDensSpect2D_',num2str(I),'_mc',num2str(JJ),'.mat'],'VarDensSpect2D')
        end
       end
 
    %%%setup signal according to simulation time
    timeSimul.interval=[timeSimul.t_init:timeSimul.dt:timeSimul.t_end];
    timeSimul.Nt=length(timeSimul.interval);

    %%%%%%%%%%%%%%%%%%%Generation Method Setup%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
            
          
        insigIhat = fft2(insigI);
        if strcmp(spatial.influx.line(I).Orientation,'Vertical')
            [WW,KK]=meshgrid(omsig,dom.ky);
        else
            [WW,KK]=meshgrid(omsig,dom.kx);
        end
        
       
        
        %         figure
        %         mesh(WW,KY,abs(insigYhat))
        %view(2)
        [~,sppeakindex]=max(max(abs(insigIhat)'));
        
        nupeak(I)      = abs(omsig(sppeakindex)); % this is peak-frequency!!!!
        
        
        %     if nupeak==0
        %     [pks,locs]   = findpeaks(abs(insighat)); % this is peak-frequency!!!!
        %     if length(locs)>1
        %     nupeak       =  abs(omsig(locs(2)));
        %     else
        %     nupeak      = abs(omsig(sppeakindex+1));
        %     end
        %     end
        %    numean      = 2*((max(omsig,0).*insighat)'*insighat)/(insighat'*insighat);%%mean period
        
        if input.wave.type(I).flag ~=2
            Tp(I)    =2*pi/nupeak(I) ;
            Hs(I)    =4*sqrt(var(reshape(insigI(:,Ind1:Ind2),[],1)));
        end
        if input.wave.type(I).flag  ==  3 || input.wave.type(I).flag  ==  5
            Sww=varDensSpect;oms=halfomsig;
            Tm01(I)=2*pi*trapz(oms,Sww)/trapz(oms,oms.*Sww);%mean period
        else
            [oms,Sww]=funSP_variance_density_spectrum1d(timesig,InfluxSignal(3:end,floor(end/2)+1));
            Tm01(I)=2*pi*trapz(oms,Sww)/trapz(oms,oms.*Sww);%mean period
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%
        ramp2dorg=funG_ramp2d(input.wave,dom,timesig,spatial.influx,Tp(I),I);
        [influx.wave(I).time,influx.wave(I).eta]=funG_adjust_tsignal_to_tsimul(timesig,insigI,timeSimul,ramp2dorg);  
        ramp2d=funG_ramp2d(input.wave,dom,influx.wave(I).time,spatial.influx,Tp(I),I);
         
        ww              = omsig;
        kx              = dom.kx;
        ky              = dom.ky;
        
        Infmodel        = model.influx.type(I);
        Infdirection    = model.influx.direction(I);
        
        if strcmp(spatial.influx.line(I).Orientation,'Vertical')
            k_half         =kx(1:floor(end/2));
        else
            k_half         =ky(1:floor(end/2));
        end
        
        
        Omega         =OmExct(k_half,meandepthinflux(I),g,omAdd);
        [~,Idnu]      =max(Omega);
        k_p(I)        = interp1(Omega(1:Idnu),k_half(1:Idnu),nupeak(I),'spline');% invOmExact(2*pi/nupeak,depth);
        lambda_p(I)   = 2*pi/k_p(I);
        Cg_p(I)=UgExct(k_p(I),meandepthinflux(I),g,[]);
        
        if lambda_p(I)/meandepthinflux(I)<2
            wavecateg(I)   ={'Deep water'};
        elseif lambda_p(I)/meandepthinflux(I)>20
            wavecateg(I)   ={'Shallow water'};
        else
            wavecateg(I)   ={'Intermediate depth'};
        end
        
        if model.nonlinear>1
            if input.wave.adjzone.check==1
            adjcoef     =input.wave.adjzone.lengthfact;
            influxPosition.x=inflXY(:,1);influxPosition.y=inflXY(:,2);          
           
            if (influxPosition.x(1)==dom.X(1) && influxPosition.x(end)==dom.X(end)) || ...
                (influxPosition.y(1)==dom.Y(1) && influxPosition.y(end)==dom.Y(end))
                ChiAdjI     = funC_cfAdj2d(dom,lambda_p(I),influxPosition,Infdirection,adjcoef);
            else
               ChiAdjI     = funC_cfNonLinAdj2d(dom,lambda_p(I),influxPosition,Infdirection,adjcoef);
            end
            
            
            ChiAdj      = ChiAdj.*ChiAdjI;
            else
            ChiAdj      =1;    
            end
          
        else
            ChiAdj      =1;
        end
 
       if abs(max(depthinfluxline)-min(depthinfluxline))<1e-3
                    [influx.gen(I)]=GenerationMethod(g,model.influx,dom,UgExct,ugAdd,meandepthinflux(I),...
                        spatial.influx,influx,ramp2d,I,timeSimul,input.wave.tapered); %derivation in space for skewgam
%             [influx.gen(I)]=GenerationMethod0(g,model.influx,dom,OmExct,UgExct,ugAdd,meandepthinflux(I),...
%                 spatial.influx,influx,ramp2d,I,timeSimul,input.wave.tapered);%%derivation in time

        else
            [InterpVg]= funOprt_GroupVelInterpolation_3p(g,dom,UgExct,OmExct,omAdd,nupeak(I),...
                depthinfluxline,bath,Proj,model.dyn);
            [influx.gen(I)]=GenerationMethod_bathy(g,model.influx,dom,UgExct,ugAdd,depthinfluxline,...
                spatial.influx,influx,ramp2d,I,InterpVg,timeSimul,input.wave.tapered);%derivation in space for skewgam
%             [influx.gen(I)]=GenerationMethod_bathy0(g,model.influx,dom,OmExct,UgExct,ugAdd,depthinfluxline,...
%                 spatial.influx,influx,ramp2d,I,InterpVg,timeSimul,input.wave.tapered);%%derivation in time for skewgam
        end
   
end
    
%%%%%%%%%%Properties influx for log file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    influx.ChiAdj     =ChiAdj;
    influx.model      =model.influx;
    influx.spatial    =spatial.influx;
    influx.input      =input.wave;
    influx.par.dt     =timesig(2)-timesig(1);
    influx.par.Nt     =length(timesig);
    influx.par.type   =input.wave.type;
    influx.par.spatial=spatial.influx;
    influx.par.k_p=k_p;   %save in struct influx;
    influx.par.lambda_p=lambda_p;
    influx.par.nu_p=nupeak;
    influx.par.T_p=Tp;
    influx.par.Tm01=Tm01;%mean period
    influx.par.Hs=Hs;
    influx.par.kh=meandepthinflux.*k_p;
    influx.par.ka=k_p.*Hs/2;
    influx.par.Cp_p=nupeak./k_p;
    influx.par.Cg_p=Cg_p;
    influx.par.lambda_per_D=lambda_p./meandepthinflux;
    influx.par.meandepth=meandepthinflux;
    influx.par.category=wavecateg;
    meandepth=mean(meandepthinflux);
end

%%%Setup for Initial Value Problem
if ~strcmpi(ivp.typename,'Zero')
    
    if strcmpi(ivp.typename,'User-defined')
        if model.phiForm==1
        [ivp.eta,ivp.phi]     =funG_IVP_data2d(ivp.userdata.data,dom);
        else
        [ivp.eta,ivp.u,ivp.v] =funG_IVP_data2duv(ivp.userdata.data,dom);    
        end
        ivp.type=6;
    else
        Ampl       =ivp.A;
        sigma      =ivp.sigma;
        if strcmpi(ivp.typename,'Gaussian (along x-axis)')
            xshift     =ivp.x0;
            InitValeta=Ampl.*exp(-(dom.X-xshift).^2./(sigma^2));
            ivp.eta=repmat(InitValeta,dom.Ny,1).*dom.cfSA;
            ivp.type=2;
        elseif strcmpi(ivp.typename,'Gaussian (along y-axis)')  
            yshift     =ivp.x0;
            InitValeta=Ampl.*exp(-(dom.Y-yshift).^2./(sigma^2));
            ivp.eta=repmat(InitValeta,1,dom.Nx).*dom.cfSA;
            ivp.type=3;
        elseif strcmpi(ivp.typename,'N-waves (along x-axis)')
            xshift     =ivp.x0;
            InitValeta=Ampl.*exp(-(dom.X-xshift).^2./(sigma^2)).*(-2*(dom.X-xshift)/(sigma^2));
            InitValeta=(Ampl/(max(InitValeta)))*InitValeta;
            ivp.eta=repmat(InitValeta,dom.Ny,1).*dom.cfSA;
            ivp.type=4;
        elseif strcmpi(ivp.typename,'N-waves (along y-axis)')
            yshift     =ivp.x0;
            InitValeta=Ampl.*exp(-(dom.Y-yshift).^2./(sigma^2)).*(-2*(dom.Y-yshift)/(sigma^2));
            InitValeta=(Ampl/(max(InitValeta)))*InitValeta;
            ivp.eta=repmat(InitValeta,1,dom.Nx).*dom.cfSA;
            ivp.type=5;
        end
        if model.phiForm==1
        ivp.phi  =zeros(dom.Ny,dom.Nx);
        else
        ivp.u  =zeros(dom.Ny,dom.Nx);  ivp.v  =zeros(dom.Ny,dom.Nx);  
        end
    end

    
    if strcmpi(input.wave.option,'No')
        meandepth     = mean(mean(-dom.bathy.profile));
        spec=funSP_spec_from_images(ivp.eta,dom.X,dom.Y,meandepth,par.g);
        Tp=spec.Tp;
        nupeak=2*pi/Tp;
        Omega   =OmExct(spec.k,meandepth,g,omAdd);
        k_p= interp1(Omega,spec.k,nupeak,'spline');
        varProf=var(reshape(ivp.eta,[],1));
        Hs=4*sqrt(varProf);%spec.Hs;
        lambda_p    = 2*pi/k_p;
        ivp.par.nupeak=nupeak;
        ivp.par.Tp=Tp;
        ivp.par.Hs=Hs;
        ivp.par.lambda_p=lambda_p;
        ivp.par.meandepth=meandepth;
        ivp.par.k_p=k_p;
        ivp.par.Cg_p=funOprt_UgExact(k_p,meandepth,g,[]);
    end
else ivp.type=1;
end

%%% Setup for boundary assimilation data
if input.bdyassim.option==1
 meandepth     = mean(mean(-dom.bathy.profile));   
[bdyassim.charupdate,ChiAdj,bdyassim.propdir]=funG_BdyAssimChar(dom,input.bdyassim);
Nt=input.bdyassim.assimdata(1,3);
Nx=input.bdyassim.assimdata(2,3);
Ny=input.bdyassim.assimdata(3,3);
bdyeta=input.bdyassim.assimdata(4+(Nt-1)*Ny:3+Nt*Ny,:);
Xdat=linspace(input.bdyassim.assimdata(2,1),input.bdyassim.assimdata(2,2),Nx);    
Ydat=linspace(input.bdyassim.assimdata(3,1),input.bdyassim.assimdata(3,2),Ny);
bdyeta=funBAssim_interpEtatoSimDom(Xdat,Ydat,bdyeta,dom).*dom.cfSA;  
spec=funSP_spec_from_images(bdyeta,dom.X,dom.Y,meandepth,par.g);
Tp=spec.Tp;
nupeak=2*pi/Tp;
Omega         =OmExct(spec.k,meandepth,g,omAdd);
k_p= interp1(Omega,spec.k,nupeak,'spline');
varProf=var(reshape(bdyeta,[],1));

Hs=4*sqrt(varProf);%spec.Hs;
lambda_p    = 2*pi/k_p;
bdyassim.par.nupeak=nupeak;
bdyassim.par.Tp=Tp;
bdyassim.par.Hs=Hs;
bdyassim.par.lambda_p=lambda_p;
bdyassim.par.meandepth=meandepth;
bdyassim.par.k_p=k_p;
bdyassim.par.Cg_p=funOprt_UgExact(k_p,meandepth,g,[]);
bdyassim.par.kh=meandepth.*k_p;
bdyassim.par.ka=k_p.*Hs/2;
bdyassim.par.Cp_p=nupeak./k_p;
if lambda_p/meandepth<2
    wavecateg   ={'Deep water'};
elseif lambda_p/meandepth>20
    wavecateg   ={'Shallow water'};
else
    wavecateg   ={'Intermediate depth'};
end
bdyassim.par.category=wavecateg;       
%[bdyassim.cfSA]=funBAssim_simulationAreaChar(dom,spatial.fbl,lambda_p);    
else
bdyassim=[];   
end
par.Hs=max(Hs); %% a threshold value in ODE

if strcmpi(input.wave.option,'No')
    %%%% For Initial Value Problem
    model.influx.direction='None';
    influx.ChiAdj =1;
    timeSimul.interval=[timeSimul.t_init:timeSimul.dt:timeSimul.t_end];
    timeSimul.Nt=length(timeSimul.interval);
end

if input.bdyassim.option==1
    influx.ChiAdj =ChiAdj;
end
%influx.ChiAdjWave=ChiAdj;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%Wall as influxing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(dom.wall.option,'yes')
[dom.wall.Gam,dom.wall.SkewGam]=funW_FluxSpatGam2D_1(g,dom,UgExct,ugAdd);
if dom.wall.ReflCoef.FlagFreqDep==1
    [dom.wall.ReflCoef.fun_Om,dom.wall.ReflCoef.fun_K]=funW_ReflCoefOmega_to_RhoK(g,dom,OmExct,omAdd,ww);
end

ChiAdjW=ones(size(dom.XX));
if model.nonlinear>1
    if strcmpi(dom.wall.option,'yes')
        if dom.wall.NInfl>0
            adjcoef     =input.wave.adjzone.lengthfact;
            FlagId=[-2 -1 1 2];
            dirflag=dom.wall.bdyInfl.dirflag;
            for ii=1:4
            if FlagId(ii)==-2
            IDdir=dirflag==-2;
            Infdirection='UniB';%'BiTB';%'UniB';
            elseif FlagId(ii)==-1
            IDdir=dirflag==-1;
            Infdirection='UniL';%'BiRL';%'UniL';
            elseif FlagId(ii)==1
            IDdir=dirflag==1;
            Infdirection='UniR';%'BiRL';%'UniR';
            elseif FlagId(ii)==2
            IDdir=dirflag==2;
            Infdirection='UniT';%'BiTB';%'UniT';
            end
            IndWallN= dom.wall.bdyInfl.index(IDdir); 
            influxPosition.x=dom.XX(IndWallN);influxPosition.y=dom.YY(IndWallN);
            ChiAdjI     = funC_cfAdj2d(dom,max(lambda_p(1:end)),influxPosition,Infdirection,adjcoef) ;
            ChiAdjW      = ChiAdjW.*ChiAdjI;
            end
        end
        influx.ChiAdj     =ChiAdj.*ChiAdjW;
    end
end
influx.ChiAdjWall =ChiAdjW;

  global ITERdtW Swall_tS tprevW
  ITERdtW=1; tprevW=0;
%   Swall_tS=[0 zeros(1,length(xy_Wall(:,1))); timesig(2) zeros(1,length(xy_Wall(:,1)))];
Swall_tS=[0 zeros(1,length(dom.wall.bdyInfl.index)); timesig(2) zeros(1,length(dom.wall.bdyInfl.index))];

%    figure(111);
%    subplot(3,1,1)
%    surf(dom.XX,dom.YY,dom.wall.Gam,'edgecolor','none');
%    colorbar;view(2); axis equal
%    xlim([min(min(dom.XX)) max(max(dom.XX))]);
%    ylim([min(min(dom.YY)) max(max(dom.YY))])
%    subplot(3,1,2)
%    surf(dom.XX,dom.YY,dom.wall.SkewGam.All,'edgecolor','none');
%    colorbar;view(2); axis equal
%    xlim([min(min(dom.XX)) max(max(dom.XX))]);
%    ylim([min(min(dom.YY)) max(max(dom.YY))])
%    subplot(3,1,3)
%    surf(dom.XX,dom.YY,ChiAdjW,'edgecolor','none');
%    colorbar;view(2); axis equal
%    xlim([min(min(dom.XX)) max(max(dom.XX))]);
%    ylim([min(min(dom.YY)) max(max(dom.YY))])
end

% if model.nonlinear>1
% figure(112);
%    surf(dom.XX,dom.YY,ChiAdj,'edgecolor','none');
%    colorbar;view(2); axis equal
%    xlim([min(min(dom.XX)) max(max(dom.XX))]);
%    ylim([min(min(dom.YY)) max(max(dom.YY))])
% 
% end


