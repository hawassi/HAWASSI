function Swall_hat=fun_wallinfluxsource(influx,Oprt,par,eta,time)

global Swall_tS ITERdtW tprevW
IndWall=closest(par.x,par.wall.position);

%   
if time>tprevW %influx.timesig_insig(ITERdtW,1)%ITERdtW*par.dt
  Swall_tS(ITERdtW,1)=time;
  if par.wall.type==1
  Swall_tS(ITERdtW,2)=par.wall.refl_Coef*eta(IndWall);
  else
      
  Rho_k_Eta=Ifft(fft(eta).*par.wall.refl_Coef_k);
  Swall_tS(ITERdtW,2)=Rho_k_Eta(IndWall);    
  end
  ITERdtW=ITERdtW+1;  
  tprevW=time;
end
%Swt=eta(IndWall);

Idt=closest(Swall_tS(:,1),time); % in case t ode < tprev
if Idt==1, Idt=2; end
Swt_skew=trapz(Swall_tS(1:Idt,1),Swall_tS(1:Idt,2));

% Sws=zeros(size(influx.timesig_insig(:,1)));
% IndtInterp=closest(influx.timesig_insig(:,1),Swall_tS(end,1));
% Sws(1:IndtInterp)=interp1(Swall_tS(:,1),Swall_tS(:,2),influx.timesig_insig(1:IndtInterp,1),'spline');
% % [Swall_tS(end,1) influx.timesig_insig(IndtInterp,1) time]

%Sw_skew=cumtrapz(influx.timesig_insig(:,1),Sws);%
% Sw_skew=Ifft(fft(Sws)./(1i.*par.wall.KOm));%
%depthW=-par.bathy(IndWall);ww=freqspace(influx.timesig_insig(:,1));
%Sw_skew=Ifft(UpExact(par.wall.KOm,depthW).*fft(cumtrapz(influx.timesig_insig(:,1),Sws)));%
%Sw_skew=Ifft(UpExact(par.wall.KOm,depthW).*fft(Sws)./(1i.*ww));%
%Sw_skew=Ifft(fft(Sws)./(1i.*par.wall.KOm));%

%Swt_skew=interp1(influx.timesig_insig(:,1),Sw_skew,time,'spline');


Swall_hat=2*Swt_skew.*par.wall.gamX_skew_hat;%Swt.*par.wall.gamX_hat;%

end