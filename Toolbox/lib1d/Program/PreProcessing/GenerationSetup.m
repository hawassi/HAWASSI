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
%%%%%%%%%%%%%%%%%  Generation Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Om              = str2func(model.dispersion);
    OmExct          = str2func('OmExact');
    UpExct          = str2func('UpExact');
    UgExct          = str2func('UgExact');
    omAdd           = model.OmFun;
    ugAdd           = model.UgFun;
    g               = par.g;
    
if input.type==2 || input.type ==3
    Tp          = input.Tp;            % [s] Peak period for JS or period for harmonic
    
    Hs          = input.Hs;            % = significant wave height [m] for JS
    ampl        = input.Hs/2;
    
    dtsig       = input.dt;%Tp/30;        % choice of accuracy input signal
    t_init      = input.t_init;
    t_end       = input.t_end;
    timesig     = [t_init:dtsig:t_end];
    Ntsig       = length(timesig);
    if  mod(Ntsig,2)~=0  %Ntsig must even for input JonsWap                 
    timesig     = [t_init:dtsig:t_end-dtsig];    
    end
    omsig       = freqspace(timesig);
    if input.ramp.check==1
    nTramp      = GUImain.ramp.val;
    ramp        = max(sign((timesig-timesig(1))-nTramp*Tp),(1-cos((timesig-timesig(1))*pi/nTramp/Tp))/2);
    endramp     = fliplr(ramp);
    end
    halfomsig   = omsig(1:floor(end/2))';
end
%% Setup generation signal
if  input.type  ==  2 %harmonic
    if input.ramp.check==1
    insig       = ampl*ramp.*(cos(timesig*2*pi/Tp)).*endramp;
    else
    insig       = ampl*(cos(timesig*2*pi/Tp));    
    end
elseif input.type ==  3 %Jonswap
    gam         = input.JS_gamma;      % steepness coeff
    length_T    = timesig(end)-timesig(1);
    Nomg        = length(halfomsig);
    influx.JS_gamma=gam;
    [JS_spec,omg,omgp,domg] = fun_JONSWAP(Nomg,length_T,Tp,Hs,gam,g,'no');
    % Make signal with ramp at start and end
    % (NOTE: random phases will change every time calling this function!)
    rng('shuffle');pause(0.1); %%shuffle random generator; 
    phase       = 2*pi*rand(1,length(omg));
    rawsighat   = sqrt(JS_spec).*exp(-1i*phase);%.*(cos(phase)-1i.*sin(phase));%
    rawsighat   = [rawsighat, fliplr(conj(rawsighat))]/dtsig*2*pi*pi;
    rawsig      = Ifft(rawsighat);
    
    
    if input.ramp.check==1
    insig       = ramp.*rawsig.*endramp;
    else
    insig       = rawsig;  
    end
    insig       = Hs*insig/(4*sqrt(var(insig)));
   
   
elseif input.type     ==  4 % This is user-input
    t_init      = input.t_init;
    t_end       = input.t_end;
    dtsig       = input.dt;
    timesig     = [t_init:dtsig:t_end];
        
    insig       = interp1(input.timesig,input.usersignal,timesig,'spline','extrap');
     
     if t_init<input.timesig(1)
        IndtI=closest(timesig,input.timesig(1));
        insig(1:IndtI-1)=0;
     end
    
    if t_end>input.timesig(end)
        IndtF=closest(timesig,input.timesig(end));
        insig(IndtF+1:end)=0;
    end

    influx.tcoarse=1;
    omsig       = freqspace(timesig);
    halfomsig   = omsig(1:floor(end/2))';
    Nomg        = length(halfomsig);
  
    if input.filter.check==1
        LFreq=input.filter.LFreq; 
        HFreq=input.filter.HFreq;
       % insig=bandpass(timesig,insig,LFreq,HFreq)';
         insig=fun_bandpass(timesig,insig,LFreq,HFreq);
    end  
end

    
%%Saving influx data for input type 2,3&4
if input.type     ~=  1
    influx.dt=timesig(2)-timesig(1);
    timesig_insig=[timesig' insig'];
    influx.timesig_insig=timesig_insig;
    tempSig=timesig_insig;
    InfluxSignal=tempSig;
    save ('-v7.3',[Proj.Dir,'InfluxSignal.mat'],'InfluxSignal')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%Generation Method Setup%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ChiAdj=1;
if input.type~=1
    %% Make primitive of influx signal and plots
   % insig       = insig - mean(insig);
    insighat = fft(insig.');
    varin=var(insig.');
    orginflShat=smooth(abs(fft(insig.'-mean(insig))).^2,5);
    varcal=trapz(omsig(1:floor(end/2)),orginflShat(1:floor(end/2)));
    orginflShat=orginflShat.*varin./varcal;
    [~,sppeakindex]=max(orginflShat);
    nupeak      = abs(omsig(sppeakindex)); % this is peak-frequency!!!!
    if nupeak==0
    [pks,locs]   = findpeaks(abs(orginflShat)); % this is peak-frequency!!!! 
    if length(locs)>1
    nupeak       =  abs(omsig(locs(2)));
    else
    nupeak      = abs(omsig(sppeakindex+1));    
    end
    end
    
    numean      = 2*((max(omsig,0).*insighat)'*insighat)/(insighat'*insighat);%mean period
    
    if input.type==2 ||input.type==3
        nupeak=2*pi/Tp;
    end
    if input.type ==4
        Tp    =2*pi/nupeak;
        Hs    =4*sqrt(var(insig,1));
        nTramp =input.ramp.length;
        if input.ramp.check==1
        ramp        = max(sign((timesig-timesig(1))-nTramp*Tp),(1-cos((timesig-timesig(1))*pi/nTramp/Tp))/2);
        endramp     = fliplr(ramp);
        insig       = ramp.*insig.*endramp;
        end
    end
    Orginsig   =insig;       %%Original influx signal
    Orginsighat=insighat;
    
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Modified input signals for point-influxing && Area-short influxing
    depthinflux     = par.depth;
    ww              = omsig;
    k               = par.k;
    Infmodel        = model.influx.type;
    Infdirection    = model.influx.direction;
    om2k            = invOmExact(ww,depthinflux);% invom_interp(k,ww,depthinflux);

    insig_skew      = cumtrapz(timesig,insig);%Ifft(UpExct(om2k,depthinflux)'.*fft(cumtrapz(timesig,insig))); 
    
    GVom            = zeros(size(k));
    if strcmp(Infmodel,'Point')
        GVom            = UgExct(om2k,depthinflux,ugAdd);%%% Group Velocity
        modinsighat     = GVom.*insighat;
        modinsig        = Ifft(modinsighat);
        insig_skewhat   = fft(insig_skew);
        modinsig_skewhat= GVom'.*insig_skewhat;
        modinsig_skew   = Ifft(modinsig_skewhat);
        insig           = modinsig;
        insig_skew      = modinsig_skew;
       
    elseif  strcmp(Infmodel,'GVAreaShort')
        dx              = par.dx;
        fact            = 1/dx;%length(k)/(Xright-Xleft); %=1/dx
        GVom            = UgExct(om2k,depthinflux,ugAdd)*fact;%%% needs Ug to be defined
        GVom            = GVom*sqrt(g*depthinflux)/GVom(1)/dx;
        %%%these spatial factor (1/dx) is used to compensate the factor in gamX
        %%%in InfluxSource function       
    end
    k_half      =k(1:floor(end/2));
    
    if strcmp(model.dispersion,'OmKdV')||strcmp(model.dispersion,'OmBBM')
    Omega       =OmExct(k_half,depthinflux,omAdd);
    [~,Idnu]    =max(Omega); 
    k_p         = interp1(Omega(1:Idnu),k_half(1:Idnu),nupeak,'spline');% invOmExact(2*pi/nupeak,depth);
    else
    k_p         =  interp1(OmExct(k_half,depthinflux,omAdd),k_half,nupeak,'spline');%interp1(OmExct(k_half,depthinflux,omAdd),k_half,nupeak,'spline');%
    end
    
    lambda_p    = 2*pi/k_p;
    adjcoef     = bath.influx_AdjZone;
    if model.nonlinear>1
        ChiAdj      = cfAdj(x,par.Xinflux,lambda_p,adjcoef,Infdirection); %nonlinear Adjustment
        if strcmpi(par.wall.presence,'yes')
        ChiAdj      = ChiAdj.*cfAdj(x,par.wall.position,lambda_p,adjcoef,par.wall.wallInfdir);
        end
    else
        ChiAdj      =1;
    end
    
    
    [INsig,gamX,INsig_skew,gamX_skew,insig,insig_skew]=InfluxSource(par,...
        timesig,insig,insig_skew,Infmodel,OmExct,UpExct,UgExct,om2k,GVom,lambda_p,omAdd,ugAdd);
    
    spatgamX    = Ifft(gamX);
end   

%%%Setup for Initial Value Problem
if IVP.type~=1
    
    if IVP.type==4
        InitVal     =IVP_data(IVP.data,x);
    else
        Ampl       =IVP.A;
        xshift     =IVP.x0;
        sigma      =IVP.lambda;
        if IVP.type==2
            InitVal.eta=Ampl.*exp(-(x'-xshift).^2./(sigma^2));
        elseif IVP.type==3
            InitVal.eta=Ampl.*exp(-(x'-xshift).^2./(sigma^2)).*(-2*(x'-xshift)/(sigma^2));
            InitVal.eta=(Ampl/(max(InitVal.eta)))*InitVal.eta;
        end
        InitVal.u  =0;
    end
    zeta0_hat    =fft(InitVal.eta.*par.cfSA');
    zu0_hat      =fft(InitVal.u.*par.cfSA');
%     figure;
%     plot(x,InitVal.u,'r',x,InitVal.eta,'b')
    IVP.zeta0_hat=zeta0_hat;
    IVP.zu0_hat  =zu0_hat;%
    
    if input.type==1
    insighat=zeta0_hat;insig=Ifft(insighat);
    depthinflux    = -par.bathy(closest(par.x,par.Xinflux));
    
    [~,sppeakindex]=max(abs(insighat));
    k           =par.k;
    k_p         = k(sppeakindex);
    if k_p<0,  k_p=k(sppeakindex+1); end
    nupeak      = OmExct(k_p,depthinflux,omAdd); % Peak-frequency
    Tp=2*pi/nupeak;
    Hs=4*sqrt(var(insig,1));
    lambda_p    = 2*pi/k_p;
    end
end
if input.type==1 && IVP.type~=1
    %%%% For Initial Value Problem
    INsig=0;gamX=0;INsig_skew=0;gamX_skew=0;insig=0;insig_skew=0;
    Infdirection ='None';ChiAdj=ones(length(x),1);
   
    if strcmpi(par.wall.presence,'yes')
        if model.nonlinear>1
            ChiAdj      = cfAdj(x,par.wall.position,lambda_p,bath.influx_AdjZone,par.wall.wallInfdir);
        end
    end
    model.influx.direction=Infdirection;
    dtsig       = input.dt;%Tp/30;        % choice of accuracy input signal
    t_init      = input.t_init;
    t_end       = 4*ceil(input.t_end/4);
    timesig     = [t_init:dtsig:t_end-dtsig] ;
    ww          = freqspace(timesig);
    om2k        = invOmExact(ww,depthinflux);% invom_interp(k,ww,depthinflux);

 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%Boundary assimilation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if input.assim.check==1
meandepth     = mean(mean(-par.bathy));   
[bdyassim.charupdate,ChiAdjAssim,bdyassim.propdir]=funG_BdyAssimChar1D(par,input.assim);
xdat=input.assim.data.eta(1,2:end);
bdyeta=input.assim.data.eta(floor(1+end/2),2:end);
bdyeta=interp1(xdat,bdyeta,par.x,'spline');
ChiAdj=ChiAdj.*ChiAdjAssim;

    [~,sppeakindex]=max(abs(fft(bdyeta)));
    k           =par.k;
    k_p         = k(sppeakindex);
    if k_p<0,  k_p=k(sppeakindex+1); end
    nupeak      = OmExct(k_p,meandepth,omAdd); % Peak-frequency
    Tp=2*pi/nupeak;
    Hs=4*sqrt(var(bdyeta,1));
    lambda_p    = 2*pi/k_p;    
    
    if input.type==1 && IVP.type==1
    %%%% For Initial Value Problem
    INsig=0;gamX=0;INsig_skew=0;gamX_skew=0;insig=0;insig_skew=0;
    Infdirection =bdyassim.propdir;
    end
    depthinflux=meandepth;
    model.influx.direction=Infdirection;
    dtsig       = input.dt;%Tp/30;        % choice of accuracy input signal
    t_init      = input.t_init;
    t_end       = 4*ceil(input.t_end/4);
    timesig     = [t_init:dtsig:t_end-dtsig] ;
    ww          = freqspace(timesig);
    om2k        = invOmExact(ww,depthinflux);% invom_interp(k,ww,depthinflux);

else
    bdyassim=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%wind generated wave%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
influx.wind=input.wind;
if input.wind.check==1
influx.wind.indexX=[closest(x,influx.wind.xinterv(1));closest(x,influx.wind.xinterv(2))];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%Parameter of influx that will be use in ODE%%%%%%%%%%%%%%%%%%%%%%%%
influx.gen.INsig      =INsig;
influx.gen.INsig_skew =INsig_skew;
influx.gen.gamX       =gamX;
influx.gen.gamX_skew  =gamX_skew;
influx.gen.timesig    =timesig;
influx.gen.nonlinAdj  =ChiAdj;
influx.input          =input;
%%%%%%%%%%Properties influx for log file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
influx.name=input.name;
influx.position=par.Xinflux;
influx.k_p=k_p;   %save in struct influx;
influx.lambda_p=lambda_p;
influx.nu_p=nupeak;
influx.Hs=Hs;
influx.Tp=Tp;
influx.kh=depthinflux*k_p;
influx.ka=k_p.*Hs./2;
influx.Cp_p=nupeak/k_p;
influx.Vg_p=UgExact(k_p,depthinflux,[]);
influx.lambda_per_H=lambda_p./depthinflux;
influx.depth=depthinflux;
maxDepth=max(-par.bathy);
influx.khdeep=maxDepth*k_p;
influx.lambda_per_hdeep=lambda_p./maxDepth;

if influx.lambda_per_H<2
    influx.category='Deep water';
elseif influx.lambda_per_H>20
    influx.category='Shallow water';
else
    influx.category='Intermediate depth';
end

if influx.lambda_per_hdeep<2
    influx.categoryMaxDepth='Deep water';
elseif influx.lambda_per_hdeep>20
    influx.categoryMaxDepth='Shallow water';
else
    influx.categoryMaxDepth='Intermediate depth';
end

%%%%%%%%%%%%Additional Setting for Wall%%%%%%%%%%%%
if strcmp(par.wall.presence, 'Yes')
    bf0         = par.wall.length/3;
    indxWL      = closest(par.x,par.wall.position);
    indxWR      = closest(par.x,par.wall.position+ par.wall.length);
    par.wall.cf      = (cf(x,x(indxWL+1),bf0)-cf(x,x(indxWR-1)-bf0,bf0)).';% cfSimulationArea
    par.wall.dampcoef=7*influx.Cp_p/bf0;
    par.wall.dxChi =zeros(size(x)).';
    par.wall.dxChi(indxWL)=-1/par.dx;
    par.wall.dxChi(indxWR)=1/par.dx;
    par.wall.indx=[indxWL indxWR];
     sig0=nupeak^2./par.g;
     Dwall=-par.bathy(indxWL);
    par.wall.Dwall=Dwall;
    par.wall.kappa = [1i.*k_p invOmEvanescent(sig0*Dwall,par.wall.Evmodes).'./Dwall];

    
     if par.wall.type==2
         %%%%%%%%%%%%Additional Setting for Wall dependence on frequency%%%%%%%%%%%%
        [Rho_k,Rho_Om,Refl_coef_Om,Refl_coef_k]=ReflCoefOmega_to_Rhok(k, par.wall.Dwall,OmExct,omAdd,ww,par);
        par.wall.refl_Coef=Refl_coef_Om;
        par.wall.refl_Coef_k=Refl_coef_k;
        par.wall.Om=ww;
        par.wall.Rho_k  =Rho_k;
        par.wall.Rho_Om =Rho_Om;
     end
    
  %%%%%for wall%%%%%%%%%%%%%
 if par.wall.method==1
  global ITERdtW Swall_tS tprevW
  ITERdtW=1; tprevW=0;
  Swall_tS=[0 0; timesig(2) 0];
  om2kwall=invOmExact(ww,par.wall.Dwall);%
  par.wall.Om=ww;
  par.wall.KOm=om2kwall;
  par.wall.UgKOm=UgExact(par.wall.KOm,par.wall.Dwall);
  par.wall.Cp2KOm=UpExact(par.wall.KOm,par.wall.Dwall).^2;
  par.wall.gamXom_hat     = interp1(k,abs(par.wall.gamX_skew_hat),om2kwall,'spline');
 end
else
 par.wall.cf =0; par.wall.dampcoef=0;  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%% Used in Log_file display%%%%%%%%%   
par.length_domain=par.Xright-par.Xleft;
par.t_init       =influx.gen.timesig(1);
par.t_end        =influx.gen.timesig(end);
par.dt           =influx.gen.timesig(2)-influx.gen.timesig(1);
par.cutfrac      =spatial.cutfrac;


