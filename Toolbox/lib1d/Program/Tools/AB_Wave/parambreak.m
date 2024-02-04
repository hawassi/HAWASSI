%% parameters of breaking
%these 4 parameters below have been declared in MainUser.m
%par.break.KBC       = GUImain.Break_param(1); % Kin.Br.Cond.: U/C in [0.7:1]
%par.break.TC        = GUImain.Break_param(2); % Terminal Condition uF=TC*uI; TC in [0.2:0.3]
%par.break.Tchar     = GUImain.Break_param(3); % Characteristic Time T*=Tchar*Tp; Tchar in [0.1:0.5]; Tp->peak period
%par.break.delb      = 1.2              ; % Mixing length Coef
%Another parameters are below

if strcmp(model.breaking.check,'Yes')
    par.break.KBC =model.breaking.KBC;
    par.break.TC  =model.breaking.TC;
    par.break.delb=model.breaking.delb;
    if IVP.type~=1
        par.break.Tstar = 5*sqrt(par.depth*g);%based on Kennedy
        Hs=max(Ifft(zeta0_hat));
    else
        par.break.Tstar = model.breaking.Tchar*influx.Tp;  %characteristic time breaking [Tp/10;Tp/2]
    end
    
    par.break.MaxEtaInit=(influx.Hs/2)*(0.9); % as reference to find peaks (Max)
    par.break.dt     =par.dt;   %for saving data each dt
    
    Indirection=model.influx.direction;
    Xinflux=par.Xinflux;
    adjcoef=bath.influx_AdjZone;
    lambda_p=influx.lambda_p;
    
    if strcmp(Indirection,'Uni+')
        par.break.xbstart=closest(x,Xinflux+0.1*adjcoef*lambda_p);
        if  exist('indShoreInit','var') %for runup
            par.break.xbend=indShoreInit;
        else
            par.break.xbend=length(x);
        end
        par.break.xbendL=[];par.break.xbstartR=[];
    elseif strcmp(Indirection,'Uni-')
        par.break.xbend=closest(x,Xinflux-0.1*adjcoef*lambda_p);
        if  exist('indShoreInit','var') %for runup
            par.break.xbstart=indShoreInit;
        else
            par.break.xbstart=1;
        end
        par.break.xbendL=[];par.break.xbstartR=[];
    elseif  strcmp(Indirection,'Bi')%'bi'
        par.break.xbstart=1;
        par.break.xbendL  =closest(x,Xinflux-0.1*adjcoef*lambda_p);
        par.break.xbstartR=closest(x,Xinflux+0.1*adjcoef*lambda_p);
        par.break.xbend  =length(x);
    elseif strcmp(Indirection,'None') %Init Value Problem
        par.break.xbstart=1;par.break.xbend=length(x);
        par.break.xbendL=[];par.break.xbstartR=[];
    end
    
%     par.break.char=par.cfSA;
%     par.break.char(par.cfSA<1)=0;
%     par.break.char(influx.gen.nonlinAdj<0.3)=0;
end