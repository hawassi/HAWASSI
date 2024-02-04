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
function dz = RHSShore180801(time,z,model,par,influx,IVP,Oprt,ProgStatbar)
%%%%%passing parameter%%%%%%%%%%%%%%
k=par.k; g=par.g; x=par.x;
INsig       =influx.gen.INsig;
INsig_skew  =influx.gen.INsig_skew;
gamX        =influx.gen.gamX;
gamX_skew   =influx.gen.gamX_skew;
Indirection =model.influx.direction;
cfSA        =par.cfSA;
ChiAdj      =influx.gen.nonlinAdj;
jProgressBar=ProgStatbar.jProgressBar;
statusbarObj=ProgStatbar.statusbarObj;
aal         =Oprt.aal;
H_minShore  =par.interp.H_minShore;
H_minDispInterp =par.interp.H_min1;
%H_minShore  =par.interp.H_min;
IntCoef     =par.interp.IntCoef;

                  
global flagWarnIVP iterProgBar ITERdt tprev
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

H=Chi.*eta-par.bathymin;
Upeta   = Ifft(Oprt.Cpeak.*z1);
Upu     = Ifft(Oprt.Cpeak.*z2);

HeavH=Heaviside(H-H_minShore);
indDI=find(HeavH==0,1,'first');
indDF=find(HeavH==0,1,'last');
ShoreChar=HeavH;  

if isempty(indDI)
indWI=1; indWF=n;  
indEndBreak=n; 
FlagOvertop=0; 
signProp=0;
elseif indDI==1  % propagation direction to the left
ShoreChar(indDI:indDF)=zeros;      
indWI=indDF+1; indWF=n;  
indEndBreak=indWI-10;
FlagOvertop=0;
signProp=-1;
elseif indDF==n
ShoreChar(indDI:indDF)=zeros;     
indWI=1; indWF=indDI-1;  
indEndBreak=indWF+10; 
FlagOvertop=0;
signProp=1;
elseif indDI~=1 && indDF~=n
%HeavH(H<50*H_min)=0;
indW=find(HeavH==1);
indD=find(HeavH<H_minShore);
ShoreChar=HeavH;  
indEndBreak=n; 
FlagOvertop=1;
signProp=0;
end

if isfield(par,'wall')
if strcmp(par.wall.presence,'Yes')
Wallchar    =par.wall.Rho';
indEndBreak =closest(x,par.wall.position);
else
Wallchar    =ones(n,1);  
end
end 
 ShoreChar=1;%% changed on 12/08/2019 do not use the characteristic function.
 indWI=1;indWF=n;

HeavChar=Wallchar.*ShoreChar;



if strcmp(model.dispersion,'OmSWE')   % Use direct model to avoid interpolation
    if model.nonlinear==2
    Htot        =(depth+Chi.*eta);
    deluH_hat   =fft(Htot.*Ifft(z2).*ShoreChar); 
    deletaH_hat = fft(Ifft(g*z1+fft(Chi.*u.^2)./2).*ShoreChar);
    dteta_hat     = -1i*k.*deluH_hat ;
    dtu_hat       = -1i*k.*deletaH_hat ;
    end
elseif model.nonlinear==2 
%     Cm         = Oprt.Cm; ConjCm=conj(Cm);
%     Cp         = Oprt.Cp; ConjCp=conj(Cp);
%     Cc         = Oprt.Cc; ConjCc=conj(Cc);
%     Cder_m     = Oprt.Cderm;
%     Cder_p     = Oprt.Cderp;
%     Cder_c     = Oprt.Cderc;
%         
%     gm=zeros(n,1);gp=gm;gc=gp;
%     if FlagOvertop==0
%     gm(indWI:indWF)=ppval(IntCoef.m,H(indWI:indWF));
%     gp(indWI:indWF)=ppval(IntCoef.p,H(indWI:indWF));
%     gc(indWI:indWF)=ppval(IntCoef.c,H(indWI:indWF));
%     else
%     gm(indW)=ppval(IntCoef.m,H(indW));
%     gp(indW)=ppval(IntCoef.p,H(indW));
%     gc(indW)=ppval(IntCoef.c,H(indW));
%     end 
%     
%     Cu_hat                = Op_C(Cm,gm,Cp,gp,Cc,gc,u,z2);
%     Cu                    = Ifft(Cu_hat);
%     Cu_ShoreChar          = Ifft(Cu_hat).*ShoreChar.*Wallchar;
%     Cu_ShoreChar_hat      = fft(Cu_ShoreChar);
%     ConjC_Cu_ShoreChar_hat= Op_C(ConjCm,gm,ConjCp,gp,ConjCc,gc,...
%     Cu_ShoreChar,Cu_ShoreChar_hat);
%     Cder_u_hat             = Op_C(Cder_m,gm,Cder_p,gp,Cder_c,gc,u,z2);
%     Cder_u                = Ifft(Cder_u_hat);
%     
%     deluH_hat   =  1/g.*ConjC_Cu_ShoreChar_hat;
%     gEta=g.*(eta+par.bathyplus);
%      
%     deletaH_hat =  fft(gEta)+ fft(Ifft(fft(Cder_u.*Cu)./g).*Chi.*Wallchar);
%     deletaH_hat = fft(Ifft(deletaH_hat).*ShoreChar);
%     
%     dteta_hat = -1i.*k.*(deluH_hat);
%     dtu_hat   = -1i.*k.*(deletaH_hat);

    C2m1         = Oprt.C2m1;
    C2p1         = Oprt.C2p1;
    C2c1         = Oprt.C2c1;
    C2m2         = Oprt.C2m2;
    C2p2         = Oprt.C2p2;
    C2c2         = Oprt.C2c2;
    
    gm1=zeros(n,1);gp1=gm1;gc1=gp1;
    gm2=zeros(n,1);gp2=gm2;gc2=gp2;
     
    gm1(indWI:indWF)=ppval(IntCoef.m1,H(indWI:indWF));
    gp1(indWI:indWF)=ppval(IntCoef.p1,H(indWI:indWF));
    gc1(indWI:indWF)=ppval(IntCoef.c1,H(indWI:indWF));
    gm2(indWI:indWF)=ppval(IntCoef.m2,H(indWI:indWF));
    gp2(indWI:indWF)=ppval(IntCoef.p2,H(indWI:indWF));
    gc2(indWI:indWF)=ppval(IntCoef.c2,H(indWI:indWF));
 
    CCu_hat= Op_C2r(C2m1,gm1,C2p1,gp1,C2c1,gc1,...
        C2m2,gm2,C2p2,gp2,C2c2,gc2,u,z2);
    CCu_ShoreChar    =Ifft(CCu_hat).*ShoreChar.*Wallchar;
    uChar           =u.*ShoreChar.*Wallchar;

    CCu_ShoreChar(H<=H_minDispInterp) = g*H(H<=H_minDispInterp).*uChar(H<=H_minDispInterp);
    CCu_ShoreChar_hat=fft(CCu_ShoreChar);
    
    kC2_u=Ifft(k.*CCu_ShoreChar_hat);
    
    deluH_hat   =  1/g.*CCu_ShoreChar_hat;
    
    gEta=g.*(eta+par.bathyplus);
    
    deletaH_hat = fft(gEta)+ fft(Ifft(fft(u.^2/2-(kC2_u./g).^2/2)).*Chi.*Wallchar);
    deletaH_hat = fft(Ifft(deletaH_hat).*ShoreChar);
    
    dteta_hat = -1i.*k.*(deluH_hat);
    dtu_hat   = -1i.*k.*(deletaH_hat);
    
elseif model.nonlinear>2 
    Cm2= Oprt.Cm.^2;  Cp2=Oprt.Cp.^2;  Cc2=Oprt.Cc.^2;
    Cm= Oprt.Cm;  Cp=Oprt.Cp;  Cc=Oprt.Cc;
    gm = IntCoef.m;   gp = IntCoef.p;  gc = IntCoef.c;
    C_u_hat         = Op_C(Cm,gm,Cp,gp,Cc,gc,u,z2);
    C_u_Char        = Ifft(C_u_hat).*HeavChar;C_u_Char_hat=fft(C_u_Char);
    Csq_u_Char_hat  = Op_C(Cm,gm,Cp,gp,Cc,gc,C_u_Char,C_u_Char_hat);
    Csq_u_Char      =Ifft(Csq_u_Char_hat);
    Csq_u_Char(depth<0)=-1*Csq_u_Char(depth<0);
    Csq_u_Char(indDI:indDF)=0;
    Csq_u_Char_hat=fft(Csq_u_Char);
    
    
    deluH1_hat      = 1/g.*Csq_u_Char_hat;
    deletaH1_hat    = g.*fft(eta.*ShoreChar);
    
    
    M0u_hat          = Op_M03p(k,g,Cm2,gm,Cp2,gp,Cc2,gc,z2);
    M0u              = Ifft(M0u_hat);M0u(depth<0)=-1*M0u(depth<0);%M0u_hat=fft(M0u);
    M0u(indDI:indDF)=0;
    etaM0uChar_hat   = fft(eta.*M0u.*HeavChar);
    M1etaM0uChar_hat = Op_M13p(k,g,Cm2,gm,Cp2,gp,Cc2,gc,etaM0uChar_hat);
    M1etaM0uChar     = Ifft(M1etaM0uChar_hat);M1etaM0uChar(depth<0)=-1*M1etaM0uChar(depth<0);
    M1etaM0uChar(indDI:indDF)=0;
    M1etaM0uChar_hat = fft(M1etaM0uChar);
    deluH2_hat       = fft(eta.*u.*HeavChar)+M1etaM0uChar_hat;
    deletaH2_hat     = fft((u.^2 - M0u.^2).*HeavChar)./2;
    
    if model.nonlinear==2
        dteta_hat   = -1i*k.*deluH1_hat -fft(Chi.*Ifft(1i.*k.*deluH2_hat));
        dtu_hat     = -1i*k.*deletaH1_hat -fft(Chi.*Ifft(1i.*k.*deletaH2_hat));
    else
        M0u_Char             = M0u.*HeavChar;
        dxu                  = Ifft(1i.*k.*z2);
        dxetaM0u_hat         = 1i.*k.*fft(eta.*M0u);dxetaMu=Ifft(dxetaM0u_hat);
        CpdxetaM0u_hat       = Op_C(Cm,gm,Cp,gp,Cc,gc,dxetaMu,dxetaM0u_hat);   
        CpdxetaM0u_Char      = Ifft(CpdxetaM0u_hat).*HeavChar;
        CpdxetaM0u_Char_hat  = fft(CpdxetaM0u_Char);
        M0dxetaM0uChar_hat     = -1/g.*1i.*k.*Op_C(Cm,gm,Cp,gp,Cc,gc,CpdxetaM0u_Char,CpdxetaM0u_Char_hat);
        M0dxetaM0uChar         = Ifft(M0dxetaM0uChar_hat);
        M0dxetaM0uChar(depth<0)=-1*M0dxetaM0uChar(depth<0); M0dxetaM0uChar(indDI:indDF)=0;
        etaM0dxetaM0uChar_hat  = fft(eta.*M0dxetaM0uChar);
        M1etaM0dxetaM0uWallChar_hat= Op_M1(k,g,Cm2,gm,Cp2,gp,etaM0dxetaM0uChar_hat);
        M1etaM0dxetaM0uWallChar =Ifft(M1etaM0dxetaM0uWallChar_hat);
        M1etaM0dxetaM0uWallChar(depth<0)=-1*M1etaM0dxetaM0uWallChar(depth<0);
        M1etaM0dxetaM0uWallChar(indDI:indDF)=0;
        M1etaM0dxetaM0uWallChar_hat=fft(M1etaM0dxetaM0uWallChar);
        eta2_dxu               = eta.^2.*dxu;   eta2_dxu_hat = fft(eta2_dxu);
        Cp_eta2_dxu_hat        = Op_C(Cm,gm,Cp,gp,Cc,gc,eta2_dxu,eta2_dxu_hat);
        Cp_eta2_dxu_Char       = Ifft(Cp_eta2_dxu_hat).*HeavChar;
        Cp_eta2_dxu_WallChar_hat   = fft(Cp_eta2_dxu_Char);
        M1eta2_dxuChar_hat     = -1/g.*1i.*k.*Op_C(Cm,gm,Cp,gp,Cc,gc,Cp_eta2_dxu_Char,Cp_eta2_dxu_WallChar_hat);
        M1eta2_dxuChar         =Ifft(M1eta2_dxuChar_hat);
        M1eta2_dxuChar(depth<0)=-1*M1eta2_dxuChar(depth<0);
        M1eta2_dxuChar(indDI:indDF)=0;
        M1eta2_dxuChar_hat     =fft(M1eta2_dxuChar);
        deluH3_hat          =  M1etaM0dxetaM0uWallChar_hat...
            -(1/2).*M1eta2_dxuChar_hat...
            -(1/2).*1i.*k.*fft(eta.^2.*M0u_Char);
        deletaH3_hat        = fft(M0u.*M0dxetaM0uChar)...
            +fft(M0u_Char.*eta.*dxu);
        
        if model.nonlinear==3
            dteta_hat   = -1i*k.*deluH1_hat -fft(Chi.*Ifft(1i.*k.*(deluH2_hat+deluH3_hat)));
            dtu_hat     = -1i*k.*deletaH1_hat -fft(Chi.*Ifft(1i.*k.*(deletaH2_hat+deletaH3_hat)));
        else
%             dx2eta_hat                  =1i.*k.*1i.*k.*z1;
%             dx2eta                      =Ifft(dx2eta_hat);
%             dxetaM0u_hat                =1i.*k.*fft(eta.*M0u);
%             M0dxetaM0u_hat              =Op_M0(k,g,Cm2,gm,Cp2,gp,dxetaM0u_hat);
%             M0dxetaM0u                  =Ifft(M0dxetaM0u_hat);
%             dxetaM0dxetaM0uChar_hat     =1i.*k.*fft(eta.*M0dxetaM0u.*HeavChar);
%             M0dxetaM0dxetaM0uChar_hat   =Op_M0(k,g,Cm2,gm,Cp2,gp,dxetaM0dxetaM0uChar_hat);
%             M0dxetaM0dxetaM0uChar       =Ifft(M0dxetaM0dxetaM0uChar_hat);
%             etaM0dxetaM0dxetaM0uChar_hat=fft(eta.*M0dxetaM0dxetaM0uChar);
%             M1etaM0dxetaM0dxetaM0uChar_hat=Op_M1(k,g,Cm2,gm,Cp2,gp,etaM0dxetaM0dxetaM0uChar_hat);
%             
%             dxeta2dxuChar_hat         =1i.*k.*fft(eta.^2.*dxu.*HeavChar);
%             M0dxeta2dxuChar_hat       =Op_M0(k,g,Cm2,gm,Cp2,gp,dxeta2dxuChar_hat);
%             M0dxeta2dxuChar           =Ifft(M0dxeta2dxuChar_hat);
%             etaM0dxeta2dxuChar_hat    =fft(eta.*M0dxeta2dxuChar);
%             M1etaM0dxeta2dxuChar_hat  =Op_M1(k,g,Cm2,gm,Cp2,gp,etaM0dxeta2dxuChar_hat);
%             
%             dxeta3dxM0uChar_hat      =1i.*k.*fft(eta.^3.*HeavChar.*Ifft(1i.*k.*M0u_hat));
%             M1dxeta3dxM0uChar_hat    =Op_M1(k,g,Cm2,gm,Cp2,gp,dxeta3dxM0uChar_hat);
%             
%             eta2dx2etaM0uChar_hat    =fft(HeavChar.*eta.^2.*dx2eta.*M0u);
%             M1eta2dx2etaM0uWallChar_hat  =Op_M(k,g,Cm2,gm,Cp2,gp,eta2dx2etaM0uChar_hat);
%             
%             deluH4_hat     = (1/3).*1i.*k.*fft(eta.^3.*dxu.*HeavChar)...
%                 +M1etaM0dxetaM0dxetaM0uChar_hat...
%                 +(1/2).*1i.*k.*fft(HeavChar.*eta.^2.*M0dxetaM0u)...
%                 +(1/2).*M1etaM0dxeta2dxuChar_hat...
%                 +(1/3).*M1dxeta3dxM0uChar_hat...
%                 +(1/2).*M1eta2dx2etaM0uWallChar_hat;
%             deletaH4_hat     = (1/2).*fft( -(eta.^2.*dxu.^2.*HeavChar)...
%                 -HeavChar.*M0dxetaM0u.^2 ...
%                 -2.*(M0u.*M0dxetaM0dxetaM0uChar)...
%                 -2.*(eta.*dxu.*M0dxetaM0u.*HeavChar)...
%                 -(M0u.*M0dxeta2dxuChar)...
%                 +(eta.^2.*Ifft(1i.*k.*M0u_hat).^2.*HeavChar)...
%                 -(eta.*dx2eta.*M0u.^2.*HeavChar))...
%                 -(1/4).*1i.*k.*1i.*k.*fft(eta.^2.*M0u.^2.*HeavChar);
%             
%             dteta_hat=  -1i*k.*deluH1_hat -fft(Chi.*Ifft(1i.*k.*(deluH2_hat+deluH3_hat+deluH4_hat)));
%             dtu_hat  =   -1i*k.*deletaH1_hat -fft(Chi.*Ifft(1i.*k.*(deletaH2_hat+deletaH3_hat+deletaH4_hat)));
       end
    end
end
if ~isempty(indDI)
    if FlagOvertop==0
        dampchar(indDI:indDF)=ones;
        Upeta(indDI:indDF)   =Oprt.Cpeak.*H(indDI:indDF);
        Upu (indDI:indDF)    =Oprt.Cpeak.*u(indDI:indDF);
    else  
        dampchar(indD)=ones;
        Upeta(indD)   =Oprt.Cpeak.*H(indD);
        Upu (indD)    =Oprt.Cpeak.*u(indD);
    end
end


%Breaking process
if strcmp(model.breaking.check,'Yes')  
    [B]=breaking_process(eta+par.bathyplus,u,depth,time,x,par.break,indEndBreak);
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
%     Cfn=par.friction.Cf;%%chezy coefficient [0.001;0.01]
%     indexf=par.friction.index;
%     Rf=zeros(n,1);
%     Cf=Cfn;%.^2.*g./(H(indexf)).^(1/3);
%     Rf(indexf)=-Cf.*(u(indexf).*abs(u(indexf))./H(indexf));
%     Rf(HeavH==0)=0;
%     Rf_hat=fft(Rf);
 
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


dteta_hat = (dteta_hat + Source -(1+signProp*1i.*k).*fft(Upeta.*dampchar)).*aal;
dtu_hat   = (dtu_hat+Rb_hat+Rf_hat+Rw_hat-(1+signProp*1i.*k).*fft(Upu.*dampchar)).*aal;

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
    if time>ITERdt*par.dt
        %if ITERdt==1, fig2=figure;end;
        plot(x,eta,'r',x,par.bathy,'--k',x,ShoreChar,'b:')
        title(['t=',num2str(time)]);
       % ylim([-0.2 0.2])
        ITERdt=ITERdt+1;
        pause(0.0001);
        
    end
end

