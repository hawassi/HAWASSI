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

function [dz] = RHSFlatWall(time,z,model,par,influx,IVP,Oprt,ProgStatbar)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global ITERdt iterProgBar
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

Cu_hat          = Cp.*z2;
Cu              = Ifft(Cu_hat);
Cu_Wall         = Cu.*Wallchar; 
Csq_u_Wall_hat  = Cp.*fft(Cu_Wall);
deluH1_hat      = 1/g.*Csq_u_Wall_hat;
deletaH1_hat    =g.*z1;

if model.nonlinear==1 
    deluHnonlin=0;deletaHnonlin=0;
else
    M0u_hat         = M0.*z2;  M0u=Ifft(M0u_hat);
    etaM0uWall_hat  = fft(eta.*M0u.*Wallchar);   
    M1etaM0uWall_hat= M1.*etaM0uWall_hat;
    deluH2_hat      = fft(u.*eta.*Wallchar)+M1etaM0uWall_hat;  
    deletaH2_hat    = fft((u.^2-M0u.^2).*Wallchar)/2;
    
    if model.nonlinear==2
        deluHnonlin=deluH2_hat;       deletaHnonlin=deletaH2_hat;
    else
%          M0u_Wall             = Ifft(-1/g.*1i.*k.*Csq_u_Wall_hat);
%          dxu                  = Ifft(1i.*k.*z2);
%         CpdxetaM0u_hat       = Cp*1i.*k.*fft(eta.*M0u);
%         CpdxetaM0u_Wall      = Ifft(CpdxetaM0u_hat).*Wallchar;
%         MdxetaMuWall_hat     =-1/g.*1i.*k.*Cp.*fft(CpdxetaM0u_Wall);
%         M0dxetaM0uWall       = Ifft(MdxetaMuWall_hat);
%         etaM0dxetaM0uWall_hat= fft(eta.*M0dxetaM0uWall);
%         M1etaM0dxetaM0uWall_hat= M1.*etaM0dxetaM0uWall_hat;
%         eta2_dxu             = eta.^2.*dxu;   eta2_dxu_hat = fft(eta2_dxu);
%         Cp_eta2_dxu_hat      = Cp.*eta2_dxu_hat;
%         Cp_eta2_dxu_Wall     = Ifft(Cp_eta2_dxu_hat).*Wallchar;
%         M1eta2_dxuWall_hat   =-1/g.*1i.*k.*Cp.*fft(Cp_eta2_dxu_Wall);
%         
%         deluH3_hat       = M1etaM0dxetaM0uWall_hat...
%                             -(1/2).*M1eta2_dxuWall_hat...
%                             -(1/2).*1i.*k.*fft(eta.^2.*M0u_Wall);
%         deletaH3_hat     = fft(M0u.*M0dxetaM0uWall)...
%                             +fft(M0u_Wall.*eta.*dxu);
        L=1i.*k.*M1; dxu = Ifft(1i.*k.*z2);
        M1etaLetaMuwall_hat  = M1.*fft(eta.*Ifft(L.*fft(eta.*M0u.*Wallchar)));
        M1walletaLetaMu_hat  = M1.*fft(Wallchar.*eta.*Ifft(L.*fft(eta.*M0u)));
        M1dxueta2wall_hat    = M1.*fft(dxu.*eta.^2.*Wallchar);
        dxeta2Muwall_hat     = 1i.*k.*fft(eta.^2.*M0u.*Wallchar);
        
        LMuetaWall_hat   =L.*fft(M0u.*eta.*Wallchar);
        LMuetaWall_Mu_hat=fft(Ifft(LMuetaWall_hat).*M0u); 
        wallMuLetaMu_hat =fft(Wallchar.*M0u.*Ifft(L.*fft(eta.*M0u)));
        
        deluH3_hat       = -0.5*(M1etaLetaMuwall_hat+M1walletaLetaMu_hat...
                            +M1dxueta2wall_hat+dxeta2Muwall_hat);
        deletaH3_hat     = 0.5*(LMuetaWall_Mu_hat+wallMuLetaMu_hat)+...
                           fft(dxu.*eta.*M0u.*Wallchar);
                        
        if model.nonlinear==3
            deluHnonlin=deluH2_hat+deluH3_hat;       deletaHnonlin=deletaH2_hat+deletaH3_hat;
        else
            dx2eta_hat               =1i.*k.*1i.*k.*z1;
            dx2eta                   =Ifft(dx2eta_hat);
            dxetaM0u_hat             =1i.*k.*fft(eta.*M0u);
            M0dxetaM0u_hat           =M0.*dxetaM0u_hat;
            M0dxetaM0u               =Ifft(M0dxetaM0u_hat);
            dxetaM0dxetaM0uWall_hat  =1i.*k.*fft(eta.*M0dxetaM0u.*Wallchar);
            M0dxetaM0dxetaM0uWall_hat=M0.*dxetaM0dxetaM0uWall_hat;
            M0dxetaM0dxetaM0uWall    =Ifft(M0dxetaM0dxetaM0uWall_hat);
            M1etaM0dxetaM0dxetaM0uWall_hat=M1.*fft(eta.*M0dxetaM0dxetaM0uWall);
            
            dxeta2dxuWall_hat        =1i.*k.*fft(eta.^2.*dxu.*Wallchar);
            M0dxeta2dxuWall_hat      =M0.*dxeta2dxuWall_hat ;
            M0dxeta2dxuWall          =Ifft(M0dxeta2dxuWall_hat);
            M1etaMdxeta2dxuWall_hat  =M1.*fft(eta.*M0dxeta2dxuWall);

            dxeta3dxM0uWall_hat      =1i.*k.*fft(eta.^3.*Wallchar.*Ifft(1i.*k.*M0u_hat));
            M1dxeta3dxM0uWall_hat    =M1.*dxeta3dxM0uWall_hat;
            
            eta2dx2etaM0uWall_hat    =fft(Wallchar.*eta.^2.*dx2eta.*M0u);
            M1eta2dx2etaM0uWall_hat  =M1.*eta2dx2etaM0uWall_hat;

          
                        
            deluH4_hat =(1/3).*1i.*k.*fft(eta.^3.*dxu.*Wallchar) ...
                +M1etaM0dxetaM0dxetaM0uWall_hat...
                +1./2.*1i.*k.*fft(Wallchar.*eta.^2.*M0dxetaM0u)...
                +1./2.*M1etaMdxeta2dxuWall_hat...
                +1./3.*M1dxeta3dxM0uWall_hat...
                +1./2.*M1eta2dx2etaM0uWall_hat;
            deletaH4_hat = -1./2.*fft(eta.^2.*dxu.^2.*Wallchar)...
                -1./2.*fft(Wallchar.*M0dxetaM0u.^2)...
                -fft(M0u.*M0dxetaM0dxetaM0uWall)...
                -fft(M0dxetaM0u.*eta.*dxu.*Wallchar)...
                -1./2.*fft(M0u.*M0dxeta2dxuWall)...
                +1./2.*fft(eta.^2.*Ifft(1i.*k.*M0u_hat).^2.*Wallchar)...
                -1./2.*fft(eta.*dx2eta.*M0u.^2.*Wallchar)...
                -1./4.*1i.*k.*1i.*k.*fft(eta.^2.*M0u.^2.*Wallchar);
            deluHnonlin   =(deluH2_hat+deluH3_hat+deluH4_hat);
            deletaHnonlin =(deletaH2_hat+deletaH3_hat+deletaH4_hat);
        end
        
    end
end
dteta_hat   = -1i*k.*deluH1_hat-fft(Ifft(1i.*k.*deluHnonlin).*Chi);
dtu_hat     = -1i.*k.*deletaH1_hat-fft(Ifft(1i.*k.*deletaHnonlin).*Chi);

if strcmp(model.breaking.check,'Yes')
    indEndBreak=closest(x,par.wall.position);
    H=(depth+eta);
    [B]=breaking_process(eta,u,depth,time,x,par.break,indEndBreak,par.t_init);
    visc=B.*par.break.delb.^2.*(H).*ifft(dteta_hat,'symmetric');
    Flux=-visc.*Ifft(dteta_hat);
    dxFlux=gradient(Flux,par.dx);
    oneperH=1./H;
    Rb_hat=fft(oneperH.*dxFlux.*Wallchar);
else
    Rb_hat=0;
end

if par.friction.check==1
H=depth+eta;
Cm=par.friction.Cf;%%manning coefficient [0.01;0.04]
indexf=par.friction.index;
Rf=zeros(n,1);
Rf(indexf)=-g*Cm.^2.*(u(indexf).*abs(u(indexf))./(H(indexf)).^(4/3));
Rf_hat=fft(Rf);
else
Rf_hat=0;    
end

if influx.wind.check==1
Rw_hat=influx.wind.coef.*z2;   
else
Rw_hat=0;
end

dteta_hat   = (dteta_hat + Source - fft(par.dampcoef.*eta.*dampchar)).*aal;
dtu_hat     = (dtu_hat   + Rb_hat +Rf_hat+Rw_hat - fft(par.dampcoef.*u.*dampchar)).*aal;

dz          = [dteta_hat; dtu_hat];


if strcmp(Visual_ODE,'On')
    if time>ITERdt*dt
        if ITERdt==1, figure;end;
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
