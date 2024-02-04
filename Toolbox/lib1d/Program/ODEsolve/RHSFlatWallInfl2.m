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
%%%%%RHS for Flat-Bottom with Wall (Frequency independent) %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dz] = RHSFlatWallInfl2(time,z,model,par,influx,IVP,Oprt,ProgStatbar)
%%%%%passing parameter%%%%%%%%%%%%%%
k=par.k; g=par.g; x=par.x;
INsig      =influx.gen.INsig;
INsig_skew =influx.gen.INsig_skew;
gamX       =influx.gen.gamX;
gamX_skew  =influx.gen.gamX_skew;
Indirection=model.influx.direction;
depth      =par.depth;
cfSA       =par.cfSA;
ChiAdj     =influx.gen.nonlinAdj;
jProgressBar=ProgStatbar.jProgressBar;
statusbarObj=ProgStatbar.statusbarObj;
aal         =Oprt.aal;
M0           =Oprt.M0;
M1           =Oprt.M1;
Wallchar    =par.wall.Rho';
WallgamX_hat  =par.wall.gamX_hat;
WallgamX_skew_hat =par.wall.gamX_skew_hat;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global ITERdt iterProgBar Swall_tS ITERdtW tprev
Visual_ODE='Off';
global flagWarnIVP
if IVP.type~=1
if flagWarnIVP==1
    w = warning('query','last');  %turn off warning odeplot
    warning('off',w.identifier);  %for initial value problem case
    flagWarnIVP=0;
end
end

if par.IP.check==1
    global dteta dtphihat iterInterior iterInteriordt timeIP uIP etaIP
end


%%Progress bar
if time>par.t_init+iterProgBar*par.dt*floor(0.1*par.Nt) %update every 10%
    set(jProgressBar,'Maximum',par.Nt, 'Value',iterProgBar*floor(0.1*par.Nt));
    jProgressBar.setStringPainted( true );
    ETA=remain_time(iterProgBar*floor(0.1*par.Nt),par.Nt);
    statusbarObj.setText(['time remaining=', num2str(ETA)]);
    iterProgBar=iterProgBar+1;
end


dampchar    = (1-cfSA)';

Chi         = ChiAdj;
Cp          = Oprt.Upd;

if influx.input.type~=1
if strcmp(Indirection, 'Uni+')
    Source       = ppval(INsig,time).*gamX - ppval(INsig_skew,time).*gamX_skew;
elseif strcmp(Indirection, 'Uni-')
    Source       = ppval(INsig,time).*gamX + ppval(INsig_skew,time).*gamX_skew;
elseif strcmp(Indirection, 'Bi')
    Source       = 2*ppval(INsig,time).*gamX ;
end
else
    Source=0;
end


n       = length(k);
z1      = z(1:n);      % etahat
z2      = z(n+1:2*n);  % uhat
eta     = Ifft(z1);
u       = Ifft(z2);



L0phi_hat         =Oprt.M0.*z2; 
L0phi=Ifft(L0phi_hat);
IndWall=closest(x,par.wall.position);
if time>influx.timesig_insig(ITERdtW,1)%ITERdtW*par.dt
  Swall_tS(ITERdtW,1)=time;
  Swall_tS(ITERdtW,2)=eta(IndWall);
  ITERdtW=ITERdtW+1;  
end
% tprev=time;

Swt=eta(IndWall);

Sws=zeros(size(influx.timesig_insig(:,1)));
IndtInterp=closest(influx.timesig_insig(:,1),Swall_tS(end,1));
Sws(1:IndtInterp)=interp1(Swall_tS(:,1),Swall_tS(:,2),influx.timesig_insig(1:IndtInterp,1),'spline');
% [Swall_tS(end,1) influx.timesig_insig(IndtInterp,1) time]

Sw_skew=cumtrapz(influx.timesig_insig(:,1),Sws);%Ifft(fft(Sws)./(1i.*par.wall.KOm));%
% assignin('base','Swall_tS',Swall_tS)
% assignin('base','Sw_skew',Sw_skew)

%Swt_skew=trapz(Swall_tS(:,1),Swall_tS(:,2));



if strcmpi(par.wall.gamX_Check,'Point')    
Sw =  Ifft(fft(Sws).*par.wall.UgKOm);
Swt=  interp1(influx.timesig_insig(:,1),Sw,time,'spline');
Sw_skew=Ifft(fft(Sw_skew).*par.wall.UgKOm);
end

Swt_skew=interp1(influx.timesig_insig(:,1),Sw_skew,time,'spline');
Swall_hat=Swt.*par.wall.gamX_hat;%2*Swt_skew.*par.wall.gamX_skew_hat;%%

delphiH1_hat     = fft(L0phi)+Swall_hat;
dxdeletaH1_hat    =  g.*1i.*k.*fft(eta);

dteta_hat   = (delphiH1_hat + Source - fft(par.dampcoef.*eta.*dampchar)).*aal;
dtphi_hat     = (-dxdeletaH1_hat- fft(par.dampcoef.*u.*dampchar)).*aal;

dz          = [dteta_hat; dtphi_hat];



if strcmp(Visual_ODE,'On')
    if time>ITERdt*par.dt
        if ITERdt==1, figure;end
        plot(x,eta,'r',x,Wallchar,'k')
        title(gca,['t=',num2str(time)]);
        ylim([-0.2 0.2])
        pause(0.0001);
        ITERdt=ITERdt+1;
    end
end

if par.IP.check==1
    if time>=par.IP.time(1) && time<=par.IP.time(2)
        if time>=par.IP.time(1)+(iterInteriordt-1)*par.dt*par.IP.time(3)
            dteta(iterInterior,:)=Ifft(dz(1:n));
            dtphihat(iterInterior,:)=dz(n+1:end)./(1i.*k);
            timeIP(iterInterior,1)=time;
            uIP(iterInterior,:)   =u;
            etaIP(iterInterior,:) =eta;
            if any(k==0)
                dtphihat(isnan(dtphihat(iterInterior,:)))=0;
            end
            iterInterior=iterInterior+1;
            iterInteriordt=iterInteriordt+1;
        end
    end
end
