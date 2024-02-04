function [Swall_hat]=funW_influxsource_2D(influx,dom,eta,time)

global Swall_tS ITERdtW tprevW
IndWall=dom.wall.bdyInfl.index;
dirflag=dom.wall.bdyInfl.dirflag;
Npw=length(IndWall);
if time>tprevW %influx.timesig_insig(ITERdtW,1)%ITERdtW*par.dt
    Swall_tS(ITERdtW,1)=time;
   
       Refl_Eta=dom.wall.ReflCoef.uniform.*eta;
       if dom.wall.ReflCoef.FlagFreqDep==1
       for JJ=1:dom.wall.N
           if dom.wall.ReflCoef.FreqDep(JJ).flag==1
           Refl_Eta=Refl_Eta+dom.wall.ReflCoef.FreqDep(JJ).char...
                    .*funC_ifft2(fft2(eta).*dom.wall.ReflCoef.fun_K(JJ).coef);
           end
       end
       end
       Swall_tS(ITERdtW,2:Npw+1)=Refl_Eta(IndWall);

    ITERdtW=ITERdtW+1;
    tprevW=time; 
end
% if strcmpi(dom.wall.spat.line(1).Orientation,'Vertical')
% Swt=repmat(eta(IndWall),1,dom.Nx);
% else
% Swt=repmat(eta(IndWall),1,dom.Ny)';
% end

Swall_hat=0;
Swt_skew=zeros(size(dom.XX));
FlagId=[-2 -1 1 2];

Idt=closest(Swall_tS(:,1),time); % in case t ode < tprev
if Idt==1, Idt=2; end
Swt_skewLine=trapz(Swall_tS(1:Idt,1),Swall_tS(1:Idt,2:end));    
    

for ii=1:4
    SkewGam=zeros(size(dom.XX));
    if FlagId(ii)==-2
      IDdir=dirflag==-2;
      IndWallN= IndWall(IDdir);  
      SigTemp = zeros(length(dom.X),1);
      IndXwall= funC_closest(dom.X,dom.XX(IndWallN));
      SigTemp(IndXwall)=Swt_skewLine(IDdir);
      Swt_skew=repmat(SigTemp,1,dom.Ny)';
      ffact   =dom.dx./(2*pi);
      SkewGam=dom.wall.SkewGam.min2;
    elseif FlagId(ii)==-1
      IDdir=dirflag==-1;
      IndWallN= IndWall(IDdir);  
      SigTemp = zeros(length(dom.Y),1);
      IndXwall= funC_closest(dom.Y,dom.YY(IndWallN));
      SigTemp(IndXwall)=Swt_skewLine(IDdir);      
      Swt_skew=repmat(SigTemp,1,dom.Nx);
      ffact   =dom.dy./(2*pi);
      SkewGam=dom.wall.SkewGam.min1;
    elseif FlagId(ii)==1
      IDdir=dirflag==1;
      IndWallN= IndWall(IDdir);  
      SigTemp = zeros(length(dom.Y),1);
      IndXwall= funC_closest(dom.Y,dom.YY(IndWallN));
      SigTemp(IndXwall)=Swt_skewLine(IDdir);
      Swt_skew=repmat(SigTemp,1,dom.Nx);
      ffact   =dom.dy./(2*pi);
      SkewGam=dom.wall.SkewGam.plus1;
    elseif FlagId(ii)==2
      IDdir=dirflag==2;
      IndWallN= IndWall(IDdir);  
      SigTemp = zeros(length(dom.X),1);
      IndXwall= funC_closest(dom.X,dom.XX(IndWallN));
      SigTemp(IndXwall)=Swt_skewLine(IDdir);
      Swt_skew=repmat(SigTemp,1,dom.Ny)';
      ffact   =dom.dx./(2*pi); 
      SkewGam=dom.wall.SkewGam.plus2;
    else
      ffact   =0;  
    end
   
    
%     assignin('base','SkewGam',SkewGam);
%     assignin('base','Swt_skew',Swt_skew);
%     assignin('base','IndWallN',IndWallN);
    
    % Sws=zeros(size(influx.timesig_insig(:,1)));
    % IndtInterp=closest(influx.timesig_insig(:,1),Swall_tS(end,1));
    % Sws(1:IndtInterp)=interp1(Swall_tS(:,1),Swall_tS(:,2),influx.timesig_insig(1:IndtInterp,1),'spline');
    % % [Swall_tS(end,1) influx.timesig_insig(IndtInterp,1) time]
    %
    % Sw_skew=cumtrapz(influx.timesig_insig(:,1),Sws);%
    % % Sw_skew=Ifft(fft(Sws)./(1i.*par.wall.KOm));%
    % depthW=-par.bathy(IndWall);ww=freqspace(influx.timesig_insig(:,1));
    % Sw_skew=Ifft(UpExact(par.wall.KOm,depthW).*fft(cumtrapz(influx.timesig_insig(:,1),Sws)));%
    % Sw_skew=Ifft(UpExact(par.wall.KOm,depthW).*fft(Sws)./(1i.*ww));%
    
    % Swt_skew=interp1(influx.timesig_insig(:,1),Sw_skew,time,'spline');
    
    Swall_hat=Swall_hat+2*fft2(Swt_skew.*SkewGam).*ffact./(2*pi);%fft2(Swt.*dom.wall.Gam)./(2*pi);%
end
end