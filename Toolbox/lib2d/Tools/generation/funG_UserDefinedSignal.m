function   [insig,timesig,omsig,halfomsig]=funG_UserDefinedSignal(inputWave,dom,spatialinflux,I)
inflX=inputWave.usersignal(I).inflX;
inflY=inputWave.usersignal(I).inflY;
Nspat=length(inflX);
Signal=inputWave.usersignal(I).eta;
Time=inputWave.usersignal(I).time;
XYinfl=spatialinflux.line(I).xy;
if Nspat==1
    if strcmpi(spatialinflux.line(I).Orientation,'Vertical')
        insigI     = repmat(Signal',dom.Ny,1)';
        Nspat      =dom.Ny;
    else
        insigI    = repmat(Signal',dom.Nx,1)';
        Nspat      =dom.Nx;
    end
else
    if strcmpi(spatialinflux.line(I).Orientation,'Vertical')
        [S,T]=meshgrid(inflY,Time);
    else
        [S,T]=meshgrid(inflX,Time);
    end
    
    SS=reshape(S,[],1);
    TT=reshape(T,[],1);
    eta=reshape(Signal,[],1);
    Fsig=scatteredInterpolant(SS,TT,eta,'nearest');
    if strcmpi(spatialinflux.line(I).Orientation,'Vertical')
        Ind1=1;%funC_closest(dom.Y,min(XYinfl(:,2)));
        Ind2=dom.Ny;%funC_closest(dom.Y,max(XYinfl(:,2)));
        
        insigI=zeros(length(Time),length(dom.Y));
        Y=dom.Y(Ind1:Ind2);
        [yy,tt]=meshgrid(Y,Time);
        tt=reshape(tt,[],1);yy=reshape(yy,[],1);
        insigItemp=Fsig(yy,tt);
        insigItemp=reshape(insigItemp,[],length(Y));
        insigI(:,Ind1:Ind2)=insigItemp;
    else
        Ind1=1;%funC_closest(dom.X,min(XYinfl(:,1)));
        Ind2=dom.Nx;%funC_closest(dom.X,max(XYinfl(:,1)));
        
        insigI=zeros(length(Time),length(dom.X));
        X=dom.X(Ind1:Ind2);
        [xx,tt]=meshgrid(X,Time);
        tt=reshape(tt,[],1);xx=reshape(xx,[],1);
        insigItemp=Fsig(xx,tt);
        insigItemp=reshape(insigItemp,[],length(X));
        insigI(:,Ind1:Ind2)=insigItemp;
    end
end


t_init      = inputWave.t_init(I);
t_end       = inputWave.t_end(I);
IndtI       = funC_closest(Time,t_init) ;
IndtF       = funC_closest(Time,t_end) ;
dtSig       = Time(2)-Time(1);
if inputWave.dt(I)>dtSig
    tcoarse     = round(inputWave.dt(I)/dtSig);
else tcoarse=1;
end
%%%additional time


if inputWave.t_end(I)-dtSig >Time(end)
    indtF       =length(Time);
    Tadd        = inputWave.t_end(I)-Time;
    Ntadd       =floor(Tadd /dtSig)-1;
    Nttot       = indtF+Ntadd;
    tempSig     =zeros(Nttot,Nspat);
    temptimeSig=tempSig;
    temptimeSig (1:indtF)   =Time(1:indtF);
    for j=1:Ntadd
        temptimeSig(indtF+j)=temptimeSig(indtF)+dtSig*j;
    end
    tempSig(1:indtF,:)      =insigI;
    inputWave.usersignal    =tempSig;
    inputWave.timesig       =temptimeSig;
    IndtF                   = funC_closest(Time,t_end) ;
end

expsig        = insigI;
insig         = expsig(IndtI:tcoarse:IndtF,:);
timesig       = Time(IndtI:tcoarse:IndtF)';
%     influx.tcoarse=tcoarse;
omsig       = funC_freqspace(timesig);
halfomsig   = omsig(1:floor(end/2));
%  Nomg        = length(halfomsig);
%     if inputWave.filter.check==1
%         LFreq=inputWave.filter.LFreq;
%         HFreq=inputWave.filter.HFreq;
%         insig=bandpass(timesig,insig,LFreq,HFreq);
%     end

