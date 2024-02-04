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
%%%%%%%%%%%%%%%%%  RHS for Varying-Bottom %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function dz = RHSBathyWallFreq(time,z,model,par,influx,IVP,Oprt,ProgStatbar)
%%%%%passing parameter%%%%%%%%%%%%%%
k=par.k; g=par.g; x=par.x;
INsig      =influx.gen.INsig;
INsig_skew =influx.gen.INsig_skew;
gamX       =influx.gen.gamX;
gamX_skew  =influx.gen.gamX_skew;
Indirection=model.influx.direction;
cfSA       =par.cfSA;
ChiAdj     =influx.gen.nonlinAdj;
jProgressBar=ProgStatbar.jProgressBar;
statusbarObj=ProgStatbar.statusbarObj;
aal         =Oprt.aal;
Rho_k       =par.wall.Rho_k;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global flagWarnIVP iterProgBar
if IVP.type~=1 %Initial value problem
    if flagWarnIVP==1
        w = warning('query','last');  %turn off warning odeplot
        warning('off',w.identifier);
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

%%%%%%%%%%%%%%%%%%%
n       = length(k);
z1      = z(1:n);      % etahat
z2      = z(n+1:2*n);  % uhat
eta     = Ifft(z1);
u       = Ifft(z2);
z1      = fft(eta);
z2      = fft(u);

indW0=par.wall.indW0;
if par.wall.SignProp==1
indWi=indW0+1;indWf=n;%closest(x,x(end)-par.fourier);
else
indWi=1;%closest(x,x(1)+par.fourier);
indWf=indW0-1;    
end


Csq_u_hat       = OpCp2_var_Wall_B(Rho_k,z2,indWi,indWf,Oprt);
deluH1_hat      = 1/g.*Csq_u_hat;
deletaH1_hat    = g.*z1;

if model.nonlinear==1
    dteta_hat   = -1i*k.*deluH1_hat;
    dtu_hat     = -1i*k.*deletaH1_hat;
else
    M0u_hat      = OpM0_var_Wall_B(k,g,Rho_k,z2,indWi,indWf,Oprt);
    M0u          = Ifft(M0u_hat);
    etaM0u_hat   = fft(eta.*M0u);
    M1etaM0u_hat = OpM1_var_Wall_B(k,g,Rho_k,etaM0u_hat,indWi,indWf,Oprt);
    deluH2_hat  = fft(eta.*u)+M1etaM0u_hat;
    deletaH2_hat= fft(u.^2 - M0u.^2)./2;
    
    
    if model.nonlinear==2
        dteta_hat   = -1i*k.*deluH1_hat -fft(Chi.*Ifft(1i.*k.*deluH2_hat));
        dtu_hat     = -1i*k.*deletaH1_hat -fft(Chi.*Ifft(1i.*k.*deletaH2_hat));
    else
        dxu                   = Ifft(1i.*k.*z2);
        dxetaM0u_hat          = 1i.*k.*fft(eta.*M0u);
        M0dxetaM0u_hat        = OpM0_var_Wall_B(k,g,Rho_k,dxetaM0u_hat,indWi,indWf,Oprt);
        M0dxetaM0u            = Ifft(M0dxetaM0u_hat);
        etaiftM0dxetaM0u_hat  = fft(eta.*M0dxetaM0u);
        M1etaiftM0dxetaM0u_hat= OpM1_var_Wall_B(k,g,Rho_k,etaiftM0dxetaM0u_hat,indWi,indWf,Oprt);
        eta2_iftdxu_hat     = fft(eta.^2.*dxu);
        M1eta2_iftdxu_hat    = OpM1_var_Wall_B(k,g,Rho_k,eta2_iftdxu_hat,indWi,indWf,Oprt);
        
        deluH3_hat      = -( M1etaiftM0dxetaM0u_hat...
            +(1/2).*M1eta2_iftdxu_hat...
            +(1/2).*1i.*k.*fft(eta.^2.*M0u));
        deletaH3_hat    = fft(M0u.*(eta.*dxu+M0dxetaM0u));
        
        if model.nonlinear==3
            dteta_hat   = -1i*k.*deluH1_hat -fft(Chi.*Ifft(1i.*k.*(deluH2_hat+deluH3_hat)));
            dtu_hat     = -1i*k.*deletaH1_hat -fft(Chi.*Ifft(1i.*k.*(deletaH2_hat+deletaH3_hat)));
        else
            dx2eta_hat                =1i.*k.*1i.*k.*z1;
            dx2eta                    =Ifft(dx2eta_hat);
            dxetaM0dxetaM0u_hat       =1i.*k.*fft(eta.*M0dxetaM0u);
            M0dxetaM0dxetaM0u_hat     =OpM0_var_Wall_B(k,g,Rho_k,dxetaM0dxetaM0u_hat,indWi,indWf,Oprt);
            etaM0dxetaM0dxetaM0u_hat  =fft(eta.*Ifft(M0dxetaM0dxetaM0u_hat));
            M1etaM0dxetaM0dxetaM0u_hat=OpM1_var_Wall_B(k,g,Rho_k,etaM0dxetaM0dxetaM0u_hat,indWi,indWf,Oprt);
            dxeta2dxu_hat             =1i.*k.*fft(eta.^2.*dxu);
            M0dxeta2dxu_hat           =OpM0_var_Wall_B(k,g,Rho_k,dxeta2dxu_hat,indWi,indWf,Oprt);
            etaM0dxeta2dxu_hat        =fft(eta.*Ifft(M0dxeta2dxu_hat));
            M1etaM0dxeta2dxu_hat      =OpM1_var_Wall_B(k,g,Rho_k,etaM0dxeta2dxu_hat,indWi,indWf,Oprt);
            dxeta3dxM0u_hat           =1i.*k.*fft(eta.^3.*Ifft(1i.*k.*M0u_hat));
            M1dxeta3dxM0u_hat         =OpM1_var_Wall_B(k,g,Rho_k,dxeta3dxM0u_hat,indWi,indWf,Oprt);
            eta2dx2etaM0u_hat         =fft(eta.^2.*dx2eta.*M0u);
            M1eta2dx2etaM0u_hat       =OpM1_var_Wall_B(k,g,Rho_k,eta2dx2etaM0u_hat,indWi,indWf,Oprt);
            
            deluH4_hat     = -(-(1/3).*1i.*k.*fft(eta.^3.*dxu)...
                -M1etaM0dxetaM0dxetaM0u_hat...
                -(1/2).*M1etaM0dxeta2dxu_hat...
                -(1/2).*1i.*k.*fft(eta.^2.*M0dxetaM0u)...
                -(1/3).*M1dxeta3dxM0u_hat...
                -(1/2).*M1eta2dx2etaM0u_hat...
                );
            deletaH4_hat     = (1/2).*fft( -(eta.^2.*dxu.^2)...
                -M0dxetaM0u.^2 ...
                -2.*(M0u.*Ifft(M0dxetaM0dxetaM0u_hat))...
                -2.*(eta.*dxu.*M0dxetaM0u)...
                -(M0u.*Ifft(M0dxeta2dxu_hat))...
                +(eta.^2.*(Ifft(1i.*k.*M0u_hat)).^2)...
                -(eta.*dx2eta.*M0u.^2))...
                -(1/4).*1i.*k.*1i.*k.*fft(eta.^2.*M0u.^2)...
                ;
            
            dteta_hat=  -1i*k.*deluH1_hat -fft(Chi.*Ifft(1i.*k.*(deluH2_hat+deluH3_hat+deluH4_hat)));
            dtu_hat  =   -1i*k.*deletaH1_hat -fft(Chi.*Ifft(1i.*k.*(deletaH2_hat+deletaH3_hat+deletaH4_hat)));
        end
    end
end


if strcmp(model.breaking.check,'Yes')
    depth      = -par.bathy';
    [B] = breaking_process(eta,u,depth,time,x,par.break,[],par.t_init);
    visc= B.*par.break.delb.^2.*(depth+eta).*Ifft(dteta_hat);
    Rb_hat  = -fft((1./(depth+eta)).*Ifft(1i.*k.* fft(visc.*Ifft(dteta_hat)) ));
else
    Rb_hat=0;
end

if par.friction.check==1
depth= -par.bathy;
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

dteta_hat   = (dteta_hat + Source - fft(Upeta.*dampchar)).*aal;%fft(par.dampcoef.*eta.*dampchar)
dtu_hat     = (dtu_hat   + Rb_hat +Rf_hat+Rw_hat - fft(Upu.*dampchar)).*aal;%fft(par.dampcoef.*u.*dampchar)
dz         = [dteta_hat; dtu_hat];


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