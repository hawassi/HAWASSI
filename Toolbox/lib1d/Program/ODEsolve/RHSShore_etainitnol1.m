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
%%%%%%%%%%%%%%%%%  RHS for RunUp  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dz = RHSShore_etainitnol1(time,z,model,par,influx,IVP,Oprt,ProgStatbar)
%%%%%passing parameter%%%%%%%%%%%%%%
k=par.k; g=par.g; x=par.x;
INsig           =influx.gen.INsig;
INsig_skew      =influx.gen.INsig_skew;
gamX            =influx.gen.gamX;
gamX_skew       =influx.gen.gamX_skew;
Indirection     =model.influx.direction;
cfSA            =par.cfSA;
ChiAdj          =influx.gen.nonlinAdj;
jProgressBar    =ProgStatbar.jProgressBar;
statusbarObj    =ProgStatbar.statusbarObj;
aal             =Oprt.aal;
H_minDispInterp =par.interp.H_min1;
H_minShore      =par.interp.H_minShore;
IntCoef         =par.interp.IntCoef;

global flagWarnIVP iterProgBar ITERdt tprev Idstop
Visual_ODE='Off';
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

dampchar    = (1-cfSA)';
depth=-par.bathy;
Chi=ChiAdj;

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


%%Progress bar
if time>par.t_init+iterProgBar*par.dt*floor(0.1*par.Nt) %update every 10%
    set(jProgressBar,'Maximum',par.Nt, 'Value',iterProgBar*floor(0.1*par.Nt));
    jProgressBar.setStringPainted( true );
    ETA=remain_time(iterProgBar*floor(0.1*par.Nt),par.Nt);
    statusbarObj.setText(['time remaining=', num2str(ETA)]);
    iterProgBar=iterProgBar+1;
end


%%%%%%%%%%%%%%%%%%%
n       = length(k);
z1      = z(1:n);      % etahat
z2      = z(n+1:2*n);  % uhat
eta     = Ifft(z1);
u       = Ifft(z2);


Upeta   = Ifft(Oprt.Cpeak.*z1);
Upu     = Ifft(Oprt.Cpeak.*z2);

H=Chi.*eta-par.bathymin;

%%%Define active domain by a characteristic function
WaveChar=ones(size(eta));
WaveChar(H<H_minShore)=zeros;
IndxZerofirst=find(WaveChar==0,1,'first'); 
IndxZerolast=find(WaveChar==0,1,'last');
WaveChar(IndxZerofirst:IndxZerolast)=0;
dampchar(WaveChar==0)=ones;



if strcmpi(par.wall.presence,'yes')
    if par.wall.method==1
        Wallchar    =ones(n,1);
    else
        Wallchar    =par.wall.Rho';
    end
    indEndBreak =closest(x,par.wall.position);
else
    Wallchar    =ones(n,1);
    
    if depth(end)==0  %%% overtopping is not covered here.
    indEndBreak =IndxZerofirst;
    else
    indEndBreak =IndxZerolast;
    end
end
ShoreChar=1;%% changed on 12/08/2019 do not use the characteristic function.
HeavChar=Wallchar.*ShoreChar;

if strcmp(model.dispersion,'OmSWE')   % Use direct model to avoid interpolation
    if model.nonlinear==2
        Htot        =(depth+Chi.*eta);
        deluH_hat   =fft(Htot.*Ifft(z2).*ShoreChar);
        deletaH_hat = fft(Ifft(g*z1+fft(Chi.*u.^2)./2).*ShoreChar);
        dteta_hat     = -1i*k.*deluH_hat ;
        dtu_hat       = -1i*k.*deletaH_hat ;
    end    
elseif model.nonlinear>=2
    C2m1 = Oprt.C2m1; C2p1= Oprt.C2p1; C2c1  = Oprt.C2c1;
    C2m2 = Oprt.C2m2; C2p2 = Oprt.C2p2; C2c2 = Oprt.C2c2;
    Cm1 = sqrt(Oprt.C2m1); Cp1= sqrt(Oprt.C2p1); Cc1  = sqrt(Oprt.C2c1);
    Cm2 = sqrt(Oprt.C2m2); Cp2 = sqrt(Oprt.C2p2); Cc2 = sqrt(Oprt.C2c2);
    
    gm1 = IntCoef.m1;   gp1 = IntCoef.p1;  gc1 = IntCoef.c1;
    gm2 = IntCoef.m2;   gp2 = IntCoef.p2;  gc2 = IntCoef.c2;
    
    C_u_hat         = Op_Cr(Cm1,gm1,Cp1,gp1,Cc1,gc1,Cm2,gm2,Cp2,gp2,Cc2,gc2,u,z2);
    C_u_Char        = Ifft(C_u_hat).*HeavChar;C_u_Char_hat=fft(C_u_Char);
    Csq_u_Char_hat  = Op_Cr(Cm1,gm1,Cp1,gp1,Cc1,gc1,Cm2,gm2,Cp2,gp2,Cc2,gc2,C_u_Char,C_u_Char_hat);
    
    deluH1_hat      = 1/g.*Csq_u_Char_hat;
    gEta=g.*(eta+par.bathyplus);
%     if strcmp(influx.categoryMaxDepth,'Deep water')
%         gEta(H<H_minShore)=g.*H(H<H_minShore);
%     end
    deletaH1_hat    = fft(gEta.*ShoreChar);
    eta(depth<H_minDispInterp)=H(depth<H_minDispInterp);
     
    C2u_hat          = Op_C2r(C2m1,gm1,C2p1,gp1,C2c1,gc1,C2m2,gm2,C2p2,gp2,C2c2,gc2,u,z2);
    M0u_hat          = -1i.*k.*C2u_hat/g;
    M0u              = Ifft(M0u_hat);
    etaM0uChar       = eta.*M0u.*HeavChar;
    etaM0uChar_hat   = fft(etaM0uChar);
    C2etaM0uChar_hat = Op_C2r(C2m1,gm1,C2p1,gp1,C2c1,gc1,C2m2,gm2,C2p2,gp2,C2c2,gc2,etaM0uChar,etaM0uChar_hat);
    M1etaM0uChar_hat = -C2etaM0uChar_hat.*1i.*k./g;
    deluH2_hat       = fft(eta.*u.*HeavChar)+M1etaM0uChar_hat;
    deletaH2_hat     = fft((u.^2 - M0u.^2).*HeavChar)./2;
    
   if model.nonlinear==2
    dteta_hat   = -1i*k.*deluH1_hat -fft(Chi.*Ifft(1i.*k.*(deluH2_hat)));
    dtu_hat     = -1i*k.*deletaH1_hat -fft(Chi.*Ifft(1i.*k.*(deletaH2_hat)));
    else
        
            M0u_Char             = M0u.*HeavChar;
            dxu                  = Ifft(1i.*k.*z2);
            dxetaM0u_hat         = 1i.*k.*fft(eta.*M0u);dxetaMu=Ifft(dxetaM0u_hat);
            CpdxetaM0u_hat       = Op_Cr(Cm1,gm1,Cp1,gp1,Cc1,gc1,Cm2,gm2,Cp2,gp2,Cc2,gc2,dxetaMu,dxetaM0u_hat);
            CpdxetaM0u_Char      = Ifft(CpdxetaM0u_hat).*HeavChar;
            CpdxetaM0u_Char_hat  = fft(CpdxetaM0u_Char);
            CpCpdxetaM0uChar_hat   = Op_Cr(Cm1,gm1,Cp1,gp1,Cc1,gc1,Cm2,gm2,Cp2,gp2,Cc2,gc2,CpdxetaM0u_Char,CpdxetaM0u_Char_hat);
            M0dxetaM0uChar_hat     = -1/g.*1i.*k.*CpCpdxetaM0uChar_hat;
            M0dxetaM0uChar         = Ifft(M0dxetaM0uChar_hat);
            etaM0dxetaM0uChar      = eta.*M0dxetaM0uChar;
            etaM0dxetaM0uChar_hat  = fft(etaM0dxetaM0uChar);
            C2etaM0dxetaM0uWallChar_hat=Op_C2r(C2m1,gm1,C2p1,gp1,C2c1,gc1,C2m2,gm2,C2p2,gp2,C2c2,gc2,etaM0dxetaM0uChar,etaM0dxetaM0uChar_hat);
            M1etaM0dxetaM0uWallChar_hat=-C2etaM0dxetaM0uWallChar_hat.*1i.*k./g;
            eta2_dxu               = eta.^2.*dxu;   eta2_dxu_hat = fft(eta2_dxu);
            Cp_eta2_dxu_hat        = Op_Cr(Cm1,gm1,Cp1,gp1,Cc1,gc1,Cm2,gm2,Cp2,gp2,Cc2,gc2,eta2_dxu,eta2_dxu_hat);
            Cp_eta2_dxu_Char       = Ifft(Cp_eta2_dxu_hat).*HeavChar;
            Cp_eta2_dxu_WallChar_hat   = fft(Cp_eta2_dxu_Char);
            CpCpeta2_dxuChar_hat   = Op_Cr(Cm1,gm1,Cp1,gp1,Cc1,gc1,Cm2,gm2,Cp2,gp2,Cc2,gc2,Cp_eta2_dxu_Char,Cp_eta2_dxu_WallChar_hat);
            M1eta2_dxuChar_hat     = -1/g.*1i.*k.*CpCpeta2_dxuChar_hat;
            deluH3_hat          =  -M1etaM0dxetaM0uWallChar_hat...
                -(1/2).*M1eta2_dxuChar_hat...
                -(1/2).*1i.*k.*fft(eta.^2.*M0u_Char);
            deletaH3_hat        = fft(M0u.*M0dxetaM0uChar)...
                +fft(M0u_Char.*eta.*dxu);
        dteta_hat   = -1i*k.*deluH1_hat -fft(Chi.*Ifft(1i.*k.*(deluH2_hat+deluH3_hat)));
        dtu_hat     = -1i*k.*deletaH1_hat -fft(Chi.*Ifft(1i.*k.*(deletaH2_hat+deletaH3_hat)));
    end
end

Upeta(WaveChar==0)=Oprt.Cpeak.*H(WaveChar==0);
Upu(WaveChar==0)  =Oprt.Cpeak.*u(WaveChar==0);



%Breaking process
if strcmp(model.breaking.check,'Yes')
    [B]=breaking_process(eta,u,depth,time,x,par.break,indEndBreak,par.t_init);
    visc=B.*par.break.delb.^2.*(H).*ifft(dteta_hat,'symmetric');
    Flux=-visc.*Ifft(dteta_hat);
    dxFlux=gradient(Flux,par.dx);
    oneperH=zeros(n,1);
    oneperH(H>H_minShore)=1./H(H>H_minShore);
    Rb_hat=fft(oneperH.*dxFlux.*Wallchar);
    %         if time>30
    %         subplot(4,1,2)
    %         Rb=Ifft(Rbhat);
    %         plot(x,Rb,'c',x,visc,'--r',x,oneperH,':b',x,u,'k',x,ShoreChar,'g--')
    %         title(['time= ', num2str(time)])
    %         xlim([20;100]);ylim([-10;20]);
    %         pause(0.001)
    %         end
    
    
else
    Rb_hat=0;
end

if par.friction.check==1
    C2m1=par.friction.Cf;%%manning coefficient [0.01;0.04]
    indexf=par.friction.index;
    Rf=zeros(n,1);
    Rf(indexf)=-g*C2m1.^2.*(u(indexf).*abs(u(indexf))./(H(indexf)).^(4/3));
    Rf(HeavH==0)=0;
    Rf_hat=fft(Rf);
else
    Rf_hat=0;
end

if influx.wind.check==1
    Rw_hat=influx.wind.coef.*z2;
else
    Rw_hat=0;
end

if strcmpi(par.wall.presence,'yes') && par.wall.method==1
    Swall_hat=fun_wallinfluxsource(influx,par,eta,time);
else
    Swall_hat=0;
end


dteta_hat = (dteta_hat + Swall_hat +Source -fft(par.dampcoef.*Upeta.*dampchar)).*aal;
dtu_hat   = (dtu_hat+Rb_hat+Rf_hat+Rw_hat-fft(par.dampcoef.*Upu.*dampchar)).*aal;

dz  = [dteta_hat ; dtu_hat];

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

tprev=time;

if strcmp(Visual_ODE,'On')
    if time>ITERdt*par.dt && time>840
        if ITERdt==1, fig2=figure;end
        subplot(5,1,1)
        plot(x,eta,'r',x,par.bathy,'--k',x,ShoreChar,'b:')
        subplot(5,1,2)
        plot(x,u,'r')
        subplot(5,1,3)
        plot(x,Ifft(deletaH_hat),'r');
        subplot(5,1,4)
        plot(x,Ifft(deluH_hat),'r');
        subplot(5,1,5)
        plot(x,Ifft(dteta_hat),'r',x,Ifft(dtu_hat),'b')
        
        title(['t=',num2str(time)]);
       % ylim([-1 1])
        ITERdt=ITERdt+1;
        pause(0.0001);
        
    end
end

