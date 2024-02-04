function   [insigI,timesig,omsig,halfomsig]=funG_UserDefinedSignalScatInterp(inputWave,dom,spatialinflux,I)
t_init      = inputWave.t_init(I);
t_end       = inputWave.t_end(I);
dtSig       = inputWave.dt(I);



Xdat=inputWave.usersignal(I).inflX;
Ydat=inputWave.usersignal(I).inflY;
Nspat=length(Xdat);
Etadat=inputWave.usersignal(I).eta;
Timedat=inputWave.usersignal(I).time;
XYinfl=spatialinflux.line(I).xy;

if t_end<=Timedat(end)
timesig     = [t_init:dtSig:t_end];
else
timesig     = [t_init:dtSig:max(Timedat)];   
end

if Nspat==1
    if strcmpi(spatialinflux.line(I).Orientation,'Vertical')
        insigI     = repmat(Etadat',dom.Ny,1)';
        Nspat      =dom.Ny;
    else
        insigI    = repmat(Etadat',dom.Nx,1)';
        Nspat      =dom.Nx;
    end
else
    if strcmpi(spatialinflux.line(I).Orientation,'Vertical')
        Ind1d=funC_closest(Ydat,min(XYinfl(:,2)));
        Ind2d=funC_closest(Ydat,max(XYinfl(:,2)));
        
        [S,T]=meshgrid(Ydat(Ind1d:Ind2d),Timedat);
    else
        Ind1d=funC_closest(Xdat,min(XYinfl(:,1)));
        Ind2d=funC_closest(Xdat,max(XYinfl(:,1)));
        [S,T]=meshgrid(Xdat(Ind1d:Ind2d),Timedat);
    end
    
    SS=reshape(S,[],1);
    TT=reshape(T,[],1);
    etadat=reshape(Etadat,[],1);
    Fsig=scatteredInterpolant(SS,TT,etadat,'nearest');
    
    if strcmpi(spatialinflux.line(I).Orientation,'Vertical')
        Ind1=funC_closest(dom.Y,min(XYinfl(:,2)));
        Ind2=funC_closest(dom.Y,max(XYinfl(:,2)));
        
        insigI=zeros(length(timesig),length(dom.Y));
        Y=dom.Y(Ind1:Ind2);
        [yy,tt]=meshgrid(Y,timesig);
        tt=reshape(tt,[],1);yy=reshape(yy,[],1);
        insigItemp=Fsig(yy,tt);
        insigItemp=reshape(insigItemp,[],length(Y));
        insigI(:,Ind1:Ind2)=insigItemp;
    else
        Ind1=funC_closest(dom.X,min(XYinfl(:,1)));
        Ind2=funC_closest(dom.X,max(XYinfl(:,1)));
        
        insigI=zeros(length(timesig),length(dom.X));
         X=dom.X(Ind1:Ind2);
        [xx,tt]=meshgrid(X,timesig);
        
        tt=reshape(tt,[],1);xx=reshape(xx,[],1);
        insigItemp=Fsig(xx,tt);
        insigItemp=reshape(insigItemp,[],length(X));
        insigI(:,Ind1:Ind2)=insigItemp;
    end
end


assignin('base','insigI',insigI)


%%%additional time


if inputWave.t_end(I)-dtSig >Timedat(end)
    indtF       =length(timesig);
    Tadd        = inputWave.t_end(I)-timesig(end);
    Ntadd       =floor(Tadd /dtSig)-1;
    Nttot       = indtF+Ntadd;
    tempSig     =zeros(Nttot,Nspat);
    temptimeSig=tempSig;
    temptimeSig (1:indtF)   =timesig(1:indtF);
    for j=1:Ntadd
        temptimeSig(indtF+j)=temptimeSig(indtF)+dtSig*j;
    end
    tempSig(1:indtF,:)      =insigI;
    insigI                  =tempSig;
    timesig                 =temptimeSig;
   % IndtF                   = funC_closest(timesig,t_end) ;
end


% expsig        = insigI;
% insig         = expsig(IndtI:tcoarse:IndtF,:);
% timesig       = Timedat(IndtI:tcoarse:IndtF)';
%     influx.tcoarse=tcoarse;
omsig       = funC_freqspace(timesig);
halfomsig   = omsig(1:floor(end/2));
%  Nomg        = length(halfomsig);
%     if inputWave.filter.check==1
%         LFreq=inputWave.filter.LFreq;
%         HFreq=inputWave.filter.HFreq;
%         insig=bandpass(timesig,insig,LFreq,HFreq);
%     end
